using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Neo4j.Driver;

namespace Protein.Engine;

/// <summary>
/// Persists FSDE session state in the <c>run-registry</c> database.
///
/// Phase 8 responsibilities:
/// <list type="bullet">
///   <item>Store and reload per-family ε₀ calibration values across sessions.</item>
///   <item>Answer the <c>where was I?</c> query: active families, last run, coverage %.</item>
///   <item>Serve the <c>replay</c> query: return the RunStep chain for a given runId.</item>
/// </list>
/// </summary>
public sealed class FsdeSessionStore : IAsyncDisposable
{
    private readonly IDriver _driver;

    public FsdeSessionStore(string uri, string username, string password)
    {
        _driver = GraphDatabase.Driver(uri, AuthTokens.Basic(username, password));
    }

    // =========================================================================
    // ε₀ calibration persistence
    // =========================================================================

    /// <summary>
    /// Saves or updates the per-family ε₀ calibration value in <c>run-registry</c>.
    /// Creates a <c>:FamilyCalibration</c> node keyed by UniProt ID.
    /// </summary>
    public async Task SaveCalibrationAsync(string uniprotId, double epsilon0)
    {
        await using var session = _driver.AsyncSession(o => o.WithDatabase("run-registry"));
        await session.ExecuteWriteAsync(async tx =>
        {
            await tx.RunAsync(@"
                MERGE (c:FamilyCalibration { uniprotId: $uniprotId })
                SET c.epsilon0      = $epsilon0,
                    c.calibratedAt  = timestamp()",
                new { uniprotId, epsilon0 });
        });
        Console.WriteLine($"[FsdeSessionStore] Saved ε₀ = {epsilon0:F4} for {uniprotId}.");
    }

    /// <summary>
    /// Loads the previously persisted ε₀ for a family.
    /// Returns <c>null</c> if no calibration has been saved yet.
    /// </summary>
    public async Task<double?> LoadCalibrationAsync(string uniprotId)
    {
        await using var session = _driver.AsyncSession(o => o.WithDatabase("run-registry"));
        return await session.ExecuteReadAsync(async tx =>
        {
            var cursor = await tx.RunAsync(@"
                MATCH (c:FamilyCalibration { uniprotId: $uniprotId })
                RETURN c.epsilon0 AS eps",
                new { uniprotId });

            if (await cursor.FetchAsync())
                return (double?)cursor.Current["eps"].As<double>();

            return null;
        });
    }

    // =========================================================================
    // "where was I?" session summary
    // =========================================================================

    /// <summary>
    /// Returns a formatted FSDE session summary containing:
    /// <list type="bullet">
    ///   <item>All families with saved ε₀ calibrations.</item>
    ///   <item>Last 5 mutation runs (id, converged, ddg, timestamp).</item>
    ///   <item>Overall convergence rate across all runs.</item>
    ///   <item>Total certificates issued.</item>
    /// </list>
    /// </summary>
    public async Task<string> GetSessionSummaryAsync()
    {
        await using var session = _driver.AsyncSession(o => o.WithDatabase("run-registry"));

        var lines = new List<string>
        {
            "==============================================================================",
            " MutCert — FSDE Session Summary  (where was I?)",
            "=============================================================================="
        };

        // --- Calibrations ---
        var calibrations = await session.ExecuteReadAsync(async tx =>
        {
            var result = new List<(string uid, double eps, long at)>();
            var cursor = await tx.RunAsync(@"
                MATCH (c:FamilyCalibration)
                RETURN c.uniprotId AS uid, c.epsilon0 AS eps, c.calibratedAt AS at
                ORDER BY c.uniprotId");
            while (await cursor.FetchAsync())
                result.Add((
                    cursor.Current["uid"].As<string>(),
                    cursor.Current["eps"].As<double>(),
                    cursor.Current["at"].As<long>()));
            return result;
        });

        lines.Add("");
        lines.Add("Persisted ε₀ calibrations:");
        if (calibrations.Count == 0)
            lines.Add("  (none yet)");
        foreach (var (uid, eps, at) in calibrations)
        {
            var dt = DateTimeOffset.FromUnixTimeMilliseconds(at).LocalDateTime;
            lines.Add($"  {uid}  ε₀ = {eps:F4}  (calibrated {dt:yyyy-MM-dd HH:mm})");
        }

        // --- Last 5 runs ---
        var recentRuns = await session.ExecuteReadAsync(async tx =>
        {
            var result = new List<(string id, string mut, bool conv, double ddg, long ts)>();
            var cursor = await tx.RunAsync(@"
                MATCH (m:MutationRun)
                RETURN m.id AS id, m.mutationId AS mut, m.converged AS conv,
                       m.finalDDG AS ddg, m.timestamp AS ts
                ORDER BY m.timestamp DESC
                LIMIT 5");
            while (await cursor.FetchAsync())
                result.Add((
                    cursor.Current["id"].As<string>(),
                    cursor.Current["mut"].As<string>(),
                    cursor.Current["conv"].As<bool>(),
                    cursor.Current["ddg"].As<double>(),
                    cursor.Current["ts"].As<long>()));
            return result;
        });

        lines.Add("");
        lines.Add("Last 5 mutation runs:");
        if (recentRuns.Count == 0)
            lines.Add("  (none yet)");
        foreach (var (id, mut, conv, ddg, ts) in recentRuns)
        {
            var dt = DateTimeOffset.FromUnixTimeMilliseconds(ts).LocalDateTime;
            string status = conv ? "CERT" : "open";
            lines.Add($"  [{status}] {mut,-40} DDG={ddg:+0.000;-0.000} kcal/mol  ({dt:HH:mm:ss})");
        }

        // --- Convergence stats ---
        var stats = await session.ExecuteReadAsync(async tx =>
        {
            var cursor = await tx.RunAsync(@"
                MATCH (m:MutationRun)
                RETURN count(m) AS total,
                       sum(CASE WHEN m.converged THEN 1 ELSE 0 END) AS certs");
            if (await cursor.FetchAsync())
                return (cursor.Current["total"].As<long>(), cursor.Current["certs"].As<long>());
            return (0L, 0L);
        });

        lines.Add("");
        double rate = stats.Item1 > 0 ? 100.0 * stats.Item2 / stats.Item1 : 0.0;
        lines.Add($"Total runs: {stats.Item1}  |  Certificates issued: {stats.Item2}  |  Convergence rate: {rate:F1}%");
        lines.Add("==============================================================================");

        return string.Join(Environment.NewLine, lines);
    }

    // =========================================================================
    // Run-trace replay
    // =========================================================================

    /// <summary>
    /// Returns the full RunStep chain for a given <paramref name="runId"/>,
    /// ordered by <c>stepIndex</c>.
    /// </summary>
    public async Task<List<RunStepRecord>> GetRunTraceAsync(string runId)
    {
        await using var session = _driver.AsyncSession(o => o.WithDatabase("run-registry"));
        return await session.ExecuteReadAsync(async tx =>
        {
            var result = new List<RunStepRecord>();
            var cursor = await tx.RunAsync(@"
                MATCH (m:MutationRun { id: $runId })-[:HAS_STEP]->(s:RunStep)
                RETURN s.stepIndex AS idx, s.nodeId AS node, s.deltaE AS delta,
                       s.hopDist AS hop, s.lo AS lo, s.hi AS hi
                ORDER BY s.stepIndex",
                new { runId });

            while (await cursor.FetchAsync())
                result.Add(new RunStepRecord(
                    cursor.Current["idx"].As<int>(),
                    cursor.Current["node"].As<string>(),
                    cursor.Current["delta"].As<double>(),
                    cursor.Current["hop"].As<int>(),
                    cursor.Current["lo"].As<double>(),
                    cursor.Current["hi"].As<double>()));
            return result;
        });
    }

    /// <summary>
    /// Prints a human-readable run-trace replay to the console.
    /// </summary>
    public async Task PrintReplayAsync(string runId)
    {
        var steps = await GetRunTraceAsync(runId);

        Console.WriteLine("==============================================================================");
        Console.WriteLine($" Run-trace replay — {runId}");
        Console.WriteLine("==============================================================================");
        Console.WriteLine($" {"Step",-5} {"Node",-30} {"Hop",-4} {"ΔE (kcal/mol)",-16} {"lo",-10} {"hi",-10} Width");
        Console.WriteLine(" " + new string('-', 80));

        foreach (var s in steps)
        {
            double width = s.Hi - s.Lo;
            Console.WriteLine(
                $" {s.StepIndex,-5} {s.NodeId,-30} {s.HopDistance,-4} {s.DeltaEKcal,-16:+0.0000;-0.0000} {s.Lo,-10:F4} {s.Hi,-10:F4} {width:F4}");
        }

        Console.WriteLine("==============================================================================");
        Console.WriteLine($" Total steps replayed: {steps.Count}");
        Console.WriteLine("==============================================================================");
    }

    // =========================================================================
    // IAsyncDisposable
    // =========================================================================

    public async ValueTask DisposeAsync()
    {
        await _driver.DisposeAsync();
    }
}

/// <summary>Flat record used when returning RunStep data from Neo4j.</summary>
public sealed record RunStepRecord(
    int StepIndex,
    string NodeId,
    double DeltaEKcal,
    int HopDistance,
    double Lo,
    double Hi);
