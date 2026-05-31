using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Neo4j.Driver;

namespace Protein.Engine;

public class RunStep
{
    public int StepIndex { get; set; }
    public string NodeId { get; set; } = "";
    public double DeltaEKcal { get; set; }
    public int HopDistance { get; set; }
    public double CumulativeDDG { get; set; }
    public double Lo { get; set; }
    public double Hi { get; set; }
}

/// <summary>
/// Lightweight summary returned by <see cref="ConvergenceSupervisor.ProcessRunAsync"/>
/// so callers (e.g. <see cref="MutationAgentPool"/>) can inspect results without
/// re-querying Neo4j.
/// </summary>
public sealed record RunSummary(
    string MutationId,
    bool Converged,
    double FinalDDG,
    double FinalLo,
    double FinalHi,
    int TotalSteps);

public class ConvergenceSupervisor
{
    private readonly IDriver _driver;
    private readonly double _epsilon0;
    
    public ConvergenceSupervisor(string uri, string username, string password, double epsilon0)
    {
        _driver = GraphDatabase.Driver(uri, AuthTokens.Basic(username, password));
        _epsilon0 = epsilon0;
    }

    /// <summary>
    /// Processes the BFS trace, applies interval narrowing, and writes
    /// <c>MutationRun + RunStep + ConvergenceCertificate</c> to <c>run-registry</c>.
    ///
    /// The optional <paramref name="onStep"/> callback is invoked after each BFS step
    /// to allow live monitoring (Phase 8 — FSDE integration).
    /// </summary>
    public async Task<RunSummary> ProcessRunAsync(
        string mutationId,
        string targetResidueId,
        HashSet<string> shellHops,
        List<(string id, double propagatedDelta, int hopDistance)> bfsTrace,
        Action<RunStep>? onStep = null,
        double experimentalDdg = double.NaN)
    {
        var runSteps = new List<RunStep>();

        // ΔΔG_MJ = raw Miyazawa-Jernigan delta at the mutation site (sum of all
        // hop-1 deltas). Used to pre-seed the propagator interval per the spec.
        double ddgMj = bfsTrace
            .Where(t => t.hopDistance == 1)
            .Sum(t => t.propagatedDelta);

        var propagator = new EnergySignalPropagator(_epsilon0, ddgMj);

        int stableSteps = 0;
        double lastWidth = double.MaxValue;
        bool certificateEmitted = false;

        var evaluatedShell = new HashSet<string>();

        for (int k = 1; k <= bfsTrace.Count; k++)
        {
            var trace = bfsTrace[k - 1];

            propagator.Apply(trace.propagatedDelta);

            if (shellHops.Contains(trace.id))
                evaluatedShell.Add(trace.id);

            // Crossover recovery resets stable-step counter
            if (propagator.CrossoverOccurred)
            {
                stableSteps = 0;
                lastWidth   = double.MaxValue;
            }

            double currentWidth = propagator.Width;

            if (currentWidth <= lastWidth)
                stableSteps++;
            else
                stableSteps = 0;

            lastWidth = currentWidth;

            var stepObj = new RunStep
            {
                StepIndex      = k,
                NodeId         = trace.id,
                DeltaEKcal     = trace.propagatedDelta,
                HopDistance    = trace.hopDistance,
                CumulativeDDG  = propagator.CumulativeDDG,
                Lo             = propagator.Lo,
                Hi             = propagator.Hi
            };
            runSteps.Add(stepObj);
            onStep?.Invoke(stepObj);

            // Certificate criteria: width ≤ 2.0 kcal/mol, ≥ 5 monotone-narrowing steps,
            // all edges in the 8Å shell evaluated at least once.
            if (!certificateEmitted
                && propagator.Width  <= 2.0
                && stableSteps       >= 5
                && evaluatedShell.Count == shellHops.Count)
            {
                certificateEmitted = true;
                await WriteToRegistryAsync(mutationId, targetResidueId, runSteps, stepObj, true, experimentalDdg);
                break;
            }
        }
        
        if (!certificateEmitted && runSteps.Count > 0)
        {
            await WriteToRegistryAsync(mutationId, targetResidueId, runSteps, runSteps.Last(), false, experimentalDdg);
        }
        else if (!certificateEmitted)
        {
            Console.WriteLine($"Warning: Zero steps traced for {mutationId}");
        }

        var final = runSteps.Count > 0 ? runSteps.Last() : new RunStep();
        return new RunSummary(
            mutationId,
            certificateEmitted,
            final.CumulativeDDG,
            final.Lo,
            final.Hi,
            runSteps.Count);
    }

    private async Task WriteToRegistryAsync(string mutationId, string targetResidueId, List<RunStep> steps, RunStep finalStep, bool converged, double experimentalDdg = double.NaN)
    {
        await using var session = _driver.AsyncSession(o => o.WithDatabase("run-registry"));
        
        var compressedSteps = steps.Select(s => new Dictionary<string, object>
        {
            {"stepIndex", s.StepIndex},
            {"nodeId", s.NodeId},
            {"deltaE", s.DeltaEKcal},
            {"hopDist", s.HopDistance},
            {"lo", s.Lo},
            {"hi", s.Hi}
        }).ToList();

        string query = @"
            MERGE (m:MutationRun { id: $runId })
            SET m.mutationId = $mutationId,
                m.targetResidueId = $targetResidueId,
                m.finalDDG = $finalDDG,
                m.finalLo = $finalLo,
                m.finalHi = $finalHi,
                m.converged = $converged,
                m.totalSteps = $totalSteps,
                m.epsilon0 = $epsilon0,
                m.experimentalDdg = $experimentalDdg,
                m.timestamp = timestamp()
            
            WITH m
            UNWIND $steps AS step
            CREATE (s:RunStep {
                stepIndex: step.stepIndex,
                nodeId: step.nodeId,
                deltaE: step.deltaE,
                hopDist: step.hopDist,
                lo: step.lo,
                hi: step.hi
            })
            CREATE (m)-[:HAS_STEP]->(s)
        ";
        
        string certQuery = @"
            MATCH (m:MutationRun { id: $runId })
            CREATE (c:ConvergenceCertificate {
                runId: $runId,
                finalDDG: $finalDDG,
                finalLo: $finalLo,
                finalHi: $finalHi,
                width: $width,
                epsilon0: $epsilon0,
                issuedAtStep: $totalSteps,
                timestamp: timestamp()
            })
            CREATE (m)-[:ISSUED_CERTIFICATE]->(c)
        ";

        try
        {
            string runId = Guid.NewGuid().ToString();
            await session.ExecuteWriteAsync(async tx =>
            {
                await tx.RunAsync(query, new
                {
                    runId = runId,
                    mutationId,
                    targetResidueId,
                    finalDDG = finalStep.CumulativeDDG,
                    finalLo  = finalStep.Lo,
                    finalHi  = finalStep.Hi,
                    converged,
                    totalSteps = steps.Count,
                    epsilon0 = _epsilon0,
                    experimentalDdg = double.IsNaN(experimentalDdg) ? (object?)null : experimentalDdg,
                    steps = compressedSteps
                });
                
                if (converged)
                {
                    await tx.RunAsync(certQuery, new
                    {
                        runId = runId,
                        finalDDG = finalStep.CumulativeDDG,
                        finalLo  = finalStep.Lo,
                        finalHi  = finalStep.Hi,
                        width    = finalStep.Hi - finalStep.Lo,
                        epsilon0 = _epsilon0,
                        totalSteps = steps.Count
                    });
                }
            });
            Console.WriteLine($"[run-registry] Saved Trace → Mut: {mutationId} | Steps: {steps.Count} | Converged: {converged} | DDG: {finalStep.CumulativeDDG:F3} kcal/mol | Width: {(finalStep.Hi - finalStep.Lo):F2} | ε₀: {_epsilon0:F4}");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"[run-registry] Write failed: {ex.Message}");
        }
    }
}
