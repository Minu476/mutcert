using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using Neo4j.Driver;

namespace Protein.Engine;

// =============================================================================
// Phase 9 — Validation
//
// Responsibilities:
//   1. Generate/load stratified 80/20 holdout split (by |ΔΔG| severity stratum)
//      and fix it in data/s2648_split.json before any validation run.
//   2. For each held-out mutation:
//        a. Load the family graph (cached — ONE load per family).
//        b. Clone the graph, apply mutation, run BFS, run ConvergenceSupervisor.
//        c. Write experimental ΔΔG onto the MutationRun node in run-registry.
//        d. Check if experimental value falls within the converged interval.
//   3. Compute MJ-direct and mean-predictor baselines.
//   4. Print reliability calibration curve table.
//   5. Generate HTML report artifact.
//   6. Return structured results for consumption by the CLI.
// =============================================================================

/// <summary>Per-mutation validation outcome.</summary>
public sealed record MutationOutcome(
    string MutationId,
    double ExperimentalDdg,
    double MjDirectDdg,           // raw hop-1 ΔΔG_MJ
    bool Converged,
    double FinalDDG,
    double FinalLo,
    double FinalHi,
    double FinalWidth,
    int Steps);

/// <summary>Per-threshold reliability metrics for one family.</summary>
public sealed record CoverageRow(
    double WidthThreshold,
    int TotalHeldOut,
    int MutCertConverged,
    int MutCertCovered,
    double MutCertCoverage,       // converged & covered / converged (or 0 if none converged)
    double MjDirectCoverage,      // |MJ_direct - exp| ≤ threshold/2 / total
    double MeanCoverage,          // |family_mean - exp| ≤ threshold/2 / total
    double RecalibCoverage);      // linear-recalibrated interval coverage (converged only)

/// <summary>Aggregated result for one protein family.</summary>
public sealed record FamilyValidationResult(
    string FamilyName,
    string UniprotId,
    int TotalMutations,
    int TrainCount,
    int ValCount,
    double Epsilon0,
    double TrainMeanDdg,
    List<MutationOutcome> Outcomes,
    List<CoverageRow> CalibrationCurve,
    double PearsonR,        // Pearson correlation: predicted FinalDDG vs experimental (converged only)
    double SpearmanRho,     // Spearman ρ: rank correlation (converged only)
    double LinearSlope,     // OLS slope fit on train (mj_direct → experimental)
    double LinearIntercept); // OLS intercept

public class ValidationPipeline
{
    private readonly string _uri;
    private readonly string _username;
    private readonly string _password;

    // Map family DB name ↔ UniProt ID (canonical order for validation)
    private static readonly (string family, string uniprot)[] Families =
    [
        ("t4-lysozyme", "P00720"),
        ("ci2",         "P01053"),
        ("barnase",     "P00648")
    ];

    private static readonly double[] WidthThresholds = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];

    // One-letter → three-letter amino-acid code (shared with MutationDatasetLoader)
    private static readonly Dictionary<char, string> OneTo3 = new()
    {
        {'A', "ALA"}, {'R', "ARG"}, {'N', "ASN"}, {'D', "ASP"},
        {'C', "CYS"}, {'Q', "GLN"}, {'E', "GLU"}, {'G', "GLY"},
        {'H', "HIS"}, {'I', "ILE"}, {'L', "LEU"}, {'K', "LYS"},
        {'M', "MET"}, {'F', "PHE"}, {'P', "PRO"}, {'S', "SER"},
        {'T', "THR"}, {'W', "TRP"}, {'Y', "TYR"}, {'V', "VAL"}
    };

    public ValidationPipeline(string uri, string username, string password)
    {
        _uri      = uri;
        _username = username;
        _password = password;
    }

    // =========================================================================
    // Public API
    // =========================================================================

    /// <summary>Validate all three families and return their results.</summary>
    public async Task<List<FamilyValidationResult>> RunAllFamiliesAsync(
        string csvPath = "data/s2648/s2648.csv",
        string splitFile = "data/s2648_split.json")
    {
        var results = new List<FamilyValidationResult>();
        foreach (var (family, uniprot) in Families)
        {
            var result = await RunFamilyAsync(family, uniprot, csvPath, splitFile);
            results.Add(result);
        }
        return results;
    }

    /// <summary>Validate a single protein family.</summary>
    public async Task RunValidationAsync(string familyName)
    {
        (string family, string uniprot) = Families
            .FirstOrDefault(f => f.family == familyName);

        if (family is null)
            throw new ArgumentException($"Unknown family name '{familyName}'. Expected one of: t4-lysozyme, ci2, barnase.");

        var result = await RunFamilyAsync(family, uniprot,
            "data/s2648/s2648.csv", "data/s2648_split.json");

        PrintFamilyResult(result);
    }

    // =========================================================================
    // Core validation for one family
    // =========================================================================

    private async Task<FamilyValidationResult> RunFamilyAsync(
        string familyName,
        string uniprotId,
        string csvPath,
        string splitFile)
    {
        Console.WriteLine();
        Console.WriteLine($"{'=',-80}");
        Console.WriteLine($" Phase 9 Validation — {familyName} ({uniprotId})");
        Console.WriteLine($"{'=',-80}");

        // --- 1. Load or generate split ----------------------------------------
        var (train, val) = await LoadOrGenerateSplitAsync(familyName, uniprotId, csvPath, splitFile);

        Console.WriteLine($" Train: {train.Count}  Validation: {val.Count}");

        if (train.Count < 2)
        {
            Console.WriteLine(" SKIP: not enough training data.");
            return new FamilyValidationResult(familyName, uniprotId,
                train.Count + val.Count, train.Count, val.Count,
                1.5, 0.0, [], [], double.NaN, double.NaN, double.NaN, double.NaN);
        }

        // --- 2. Calibrate ε₀ from training set ONLY (no test-set leak) ---------
        double trainMean = train.Average(m => m.ExperimentalDdg);
        double eps0 = EpsilonCalibrator.CalibrateFromValues(
            train.Select(m => m.ExperimentalDdg), uniprotId);
        Console.WriteLine($" ε₀ = {eps0:F4}  train_mean ΔΔG = {trainMean:F3} kcal/mol");

        // --- 3. Load family graph ONCE and freeze it --------------------------
        Console.WriteLine($" Loading graph for {uniprotId}...");
        InMemoryGraph frozenGraph;
        try
        {
            frozenGraph = await GraphLoader.LoadGraphAsync(uniprotId, _uri, _username, _password);
            Console.WriteLine($" Graph: {frozenGraph.Nodes.Count} residues, {frozenGraph.Edges.Count} edges.");
        }
        catch (Exception ex)
        {
            Console.WriteLine($" ERROR loading graph: {ex.Message}");
            Console.WriteLine(" Continuing with empty graph (all mutations will produce zero-step traces).");
            frozenGraph = new InMemoryGraph();
        }

        // --- 3.6. Wild-type identity audit (Opus 4.8, 2026-06-01) ---------------
        // Each MutRecord carries the wildtype identity parsed from the mutation code.
        // Check it against the actual residue in the graph to detect silent numbering
        // mismatches (e.g. AlphaFold UniProt numbering vs. mature-protein PDB numbering).
        {
            int trainMismatch = 0, trainMissing = 0;
            foreach (var m2 in train)
            {
                if (!frozenGraph.Nodes.TryGetValue(m2.ResidueId, out var node)) { trainMissing++; continue; }
                if (!string.IsNullOrEmpty(m2.WildType3) && node.ResidueName != m2.WildType3) trainMismatch++;
            }
            int valMismatch = 0, valMissing = 0;
            var mismatches = new List<string>();
            foreach (var m2 in val)
            {
                if (!frozenGraph.Nodes.TryGetValue(m2.ResidueId, out var node)) { valMissing++; continue; }
                if (!string.IsNullOrEmpty(m2.WildType3) && node.ResidueName != m2.WildType3)
                {
                    valMismatch++;
                    if (mismatches.Count < 3)
                        mismatches.Add($"{m2.ResidueId}: expected {m2.WildType3}, got {node.ResidueName}");
                }
            }
            double trainPct = train.Count > 0 ? 100.0 * trainMismatch / train.Count : 0;
            double valPct   = val.Count   > 0 ? 100.0 * valMismatch   / val.Count   : 0;
            Console.WriteLine($" WT-audit: train {trainMismatch}/{train.Count} mismatch ({trainPct:F0}%), " +
                              $"{trainMissing} missing.");
            Console.WriteLine($"           val   {valMismatch}/{val.Count} mismatch ({valPct:F0}%), " +
                              $"{valMissing} missing.");
            foreach (var ex in mismatches) Console.WriteLine($"           e.g. {ex}");
            if (valMismatch > val.Count / 4)
                Console.WriteLine($" *** WARNING: >{valMismatch * 100 / val.Count}% val mismatches — results for {uniprotId} may be structurally misaligned. ***");
        }

        // --- 3.6. WT gate: drop any val (and train) records whose graph residue
        // does not match the stated wildtype. These are residual numbering mismatches
        // that the SeqOffset dictionary couldn't resolve (e.g. boundary edge cases).
        // Running them would inject noise: the agent would mutate the WRONG residue.
        var trainGated = train.Where(m2 =>
        {
            if (string.IsNullOrEmpty(m2.WildType3)) return true;
            return frozenGraph.Nodes.TryGetValue(m2.ResidueId, out var n2) && n2.ResidueName == m2.WildType3;
        }).ToList();
        var valGated = val.Where(m2 =>
        {
            if (string.IsNullOrEmpty(m2.WildType3)) return true;
            return frozenGraph.Nodes.TryGetValue(m2.ResidueId, out var n2) && n2.ResidueName == m2.WildType3;
        }).ToList();
        int trainDropped = train.Count - trainGated.Count;
        int valDropped   = val.Count   - valGated.Count;
        if (trainDropped > 0 || valDropped > 0)
            Console.WriteLine($" WT gate: dropped {trainDropped} train + {valDropped} val mismatches from run loop.");

        // --- 3.7. Linear recalibration: fit OLS on GATED train MJ-direct predictions ---
        // Phase 9 showed MJ energy is sign-anti-correlated with experiment.
        // Fit: experimental ≈ slope × MjDirect + intercept on the TRAIN set only,
        // then use the corrected center when checking validation coverage.
        Console.WriteLine($" Fitting linear recalibration on {trainGated.Count} training mutations (with ref-state correction)...");
        var trainMjPairs = GetTrainMjDirects(trainGated, frozenGraph);
        var (linSlope, linIntercept) = trainMjPairs.Count >= 2
            ? FitOLS(trainMjPairs)
            : (1.0, 0.0);
        string interceptSign = linIntercept >= 0 ? "+" : "";
        Console.WriteLine($" Linear model: corrected = {linSlope:+0.000;-0.000} × MJ_direct {interceptSign}{linIntercept:F3}");
        Console.WriteLine($"   (n = {trainMjPairs.Count} train pairs fitted)");

        // --- 4. Run each held-out mutation (WT-gated) -------------------------
        var outcomes = new List<MutationOutcome>();

        Console.WriteLine($" Running {valGated.Count} validation mutations (WT-gated, {valDropped} skipped)...");
        int done = 0;
        foreach (var mut in valGated)
        {
            var outcome = await RunOneMutationAsync(mut, frozenGraph, eps0);
            outcomes.Add(outcome);
            done++;
            if (done % 10 == 0 || done == valGated.Count)
                Console.Write($"\r  Progress: {done}/{valGated.Count}  ({done * 100.0 / valGated.Count:F0}%)   ");
        }
        Console.WriteLine();

        // --- 5. Build calibration curve + correlations -----------------------
        var curve = BuildCalibrationCurve(outcomes, trainMean, valGated.Count, linSlope, linIntercept);
        var convergedOutcomes = outcomes.Where(o => o.Converged).ToList();
        double pearson  = PearsonR(
            convergedOutcomes.Select(o => o.FinalDDG),
            convergedOutcomes.Select(o => o.ExperimentalDdg));
        double spearman = SpearmanRho(
            convergedOutcomes.Select(o => o.FinalDDG),
            convergedOutcomes.Select(o => o.ExperimentalDdg));

        Console.WriteLine($" Pearson r  = {pearson:+0.000;-0.000;0.000}  " +
                          $"Spearman ρ = {spearman:+0.000;-0.000;0.000}  " +
                          $"(n = {convergedOutcomes.Count} converged)");

        return new FamilyValidationResult(
            familyName, uniprotId,
            trainGated.Count + valGated.Count, trainGated.Count, valGated.Count,
            eps0, trainMean, outcomes, curve, pearson, spearman,
            linSlope, linIntercept);
    }

    // =========================================================================
    // Single mutation execution
    // =========================================================================

    private async Task<MutationOutcome> RunOneMutationAsync(
        MutRecord mut,
        InMemoryGraph frozenGraph,
        double eps0)
    {
        // Clone the frozen graph so this mutation doesn't affect others
        var graphCopy = frozenGraph.Clone();
        var agent     = new MutationAgent(graphCopy);

        List<(string id, double propagatedDelta, int hopDistance)> trace;
        HashSet<string> shellHops;

        try
        {
            (trace, shellHops) = agent.ApplyMutationAndGetTrace(mut.ResidueId, mut.MutantType3);
        }
        catch (Exception ex)
        {
            Console.Write($"\r  WARN: {mut.MutationId} trace failed: {ex.Message.Split('\n')[0]}  ");
            return new MutationOutcome(mut.MutationId, mut.ExperimentalDdg,
                0, false, 0, 0, 0, double.MaxValue, 0);
        }

        // MJ-direct baseline: sum of first-hop deltas only (no convergence),
        // then apply per-residue MJ reference-state correction (Opus 4.8):
        // subtract n_hop1 × (ref(mutant) − ref(wildtype)) to account for the
        // unfolded-state contact energy, bringing σ_MJ down to experimental scale.
        double mjDirect = trace.Where(t => t.hopDistance == 1).Sum(t => t.propagatedDelta);
        if (!string.IsNullOrEmpty(mut.WildType3))
        {
            int nHop1 = trace.Count(t => t.hopDistance == 1);
            mjDirect -= nHop1 * (MiyazawaJernigan.GetReferenceEnergy(mut.MutantType3)
                                - MiyazawaJernigan.GetReferenceEnergy(mut.WildType3));
        }

        if (trace.Count == 0)
        {
            return new MutationOutcome(mut.MutationId, mut.ExperimentalDdg,
                mjDirect, false, mjDirect, mjDirect, mjDirect, 0, 0);
        }

        var supervisor = new ConvergenceSupervisor(_uri, _username, _password, eps0);
        RunSummary summary;
        try
        {
            summary = await supervisor.ProcessRunAsync(
                mut.MutationId, mut.ResidueId, shellHops, trace,
                experimentalDdg: mut.ExperimentalDdg);
        }
        catch (Exception ex)
        {
            // Neo4j write may fail; still record what we computed locally
            Console.Write($"\r  WARN: {mut.MutationId} registry write failed: {ex.Message.Split('\n')[0]}  ");

            // Compute locally from trace
            var localProp = new EnergySignalPropagator(eps0, mjDirect);
            foreach (var t in trace) localProp.Apply(t.propagatedDelta);
            double lo  = localProp.Lo;
            double hi  = localProp.Hi;
            double w   = hi - lo;
            bool conv  = w <= 2.0 && trace.Count >= 5 && shellHops.Count > 0;
            summary    = new RunSummary(mut.MutationId, conv, localProp.CumulativeDDG, lo, hi, trace.Count);
        }

        return new MutationOutcome(
            mut.MutationId,
            mut.ExperimentalDdg,
            mjDirect,
            summary.Converged,
            summary.FinalDDG,
            summary.FinalLo,
            summary.FinalHi,
            summary.FinalHi - summary.FinalLo,
            summary.TotalSteps);
    }

    // =========================================================================
    // Calibration curve builder
    // =========================================================================

    private static List<CoverageRow> BuildCalibrationCurve(
        List<MutationOutcome> outcomes,
        double trainMean,
        int totalVal,
        double linSlope,
        double linIntercept)
    {
        var rows = new List<CoverageRow>();
        foreach (double w in WidthThresholds)
        {
            var converged = outcomes.Where(o => o.Converged && o.FinalWidth <= w).ToList();
            int covered   = converged.Count(o => o.ExperimentalDdg >= o.FinalLo && o.ExperimentalDdg <= o.FinalHi);

            double mutCertCov = converged.Count > 0 ? (double)covered / converged.Count : 0.0;

            // MJ-direct: |MjDirectDdg - experimental| ≤ w/2
            int mjCovered  = outcomes.Count(o => Math.Abs(o.MjDirectDdg - o.ExperimentalDdg) <= w / 2.0);
            double mjCov   = totalVal > 0 ? (double)mjCovered / totalVal : 0.0;

            // Mean predictor: |trainMean - experimental| ≤ w/2
            int meanCovered = outcomes.Count(o => Math.Abs(trainMean - o.ExperimentalDdg) <= w / 2.0);
            double meanCov  = totalVal > 0 ? (double)meanCovered / totalVal : 0.0;

            // Linear recalibration: shift converged interval center to
            // correctedCenter = slope * MjDirect + intercept; keep same half-width.
            // Addresses the anti-correlation by sign-correcting the MJ predictor.
            int recalibCovered = converged.Count(o =>
            {
                double correctedCenter = linSlope * o.MjDirectDdg + linIntercept;
                double halfWidth = o.FinalWidth / 2.0;
                return o.ExperimentalDdg >= correctedCenter - halfWidth
                    && o.ExperimentalDdg <= correctedCenter + halfWidth;
            });
            double recalibCov = converged.Count > 0 ? (double)recalibCovered / converged.Count : 0.0;

            rows.Add(new CoverageRow(w, totalVal, converged.Count, covered, mutCertCov, mjCov, meanCov, recalibCov));
        }
        return rows;
    }

    // =========================================================================
    // Console output
    // =========================================================================

    public static void PrintFamilyResult(FamilyValidationResult r)
    {
        Console.WriteLine();
        Console.WriteLine($"╔{'═',92}╗");
        Console.WriteLine($"║  {$"Reliability Calibration Curve — {r.FamilyName} ({r.UniprotId})",-90}║");
        Console.WriteLine($"║  {$"Train: {r.TrainCount}   Val: {r.ValCount}   ε₀ = {r.Epsilon0:F4}   train_mean = {r.TrainMeanDdg:+0.000;-0.000} kcal/mol",-90}║");
        if (!double.IsNaN(r.LinearSlope))
        {
            string iSign = r.LinearIntercept >= 0 ? "+" : "";
            Console.WriteLine($"║  {$"Recalib model: corrected = {r.LinearSlope:+0.000;-0.000} × MJ_direct {iSign}{r.LinearIntercept:F3} kcal/mol",-90}║");
        }
        Console.WriteLine($"╠{'═',92}╣");
        Console.WriteLine($"║  {"Width",-6} {"MutCert conv",-14} {"MutCert cov%",-14} {"Recalib cov%",-14} {"MJ-direct%",-12} {"Mean%",-10} {"Target",-7}║");
        Console.WriteLine($"╠{'═',92}╣");

        foreach (var row in r.CalibrationCurve)
        {
            string target = row.WidthThreshold == 2.0 ? (row.MutCertCoverage >= 0.80 ? "✅ PASS" : "❌ FAIL") : "";
            Console.WriteLine(
                $"║  {row.WidthThreshold,-6:F1} " +
                $"{$"{row.MutCertConverged}/{row.TotalHeldOut}",-14} " +
                $"{row.MutCertCoverage * 100,-14:F1}% " +
                $"{row.RecalibCoverage * 100,-14:F1}% " +
                $"{row.MjDirectCoverage * 100,-12:F1}% " +
                $"{row.MeanCoverage * 100,-10:F1}% " +
                $"{target,-7}║");
        }
        Console.WriteLine($"╚{'═',92}╝");
        if (!double.IsNaN(r.PearsonR))
            Console.WriteLine($"  Pearson r = {r.PearsonR:+0.000;-0.000;0.000}  " +
                              $"Spearman ρ = {r.SpearmanRho:+0.000;-0.000;0.000}");
    }

    // =========================================================================
    // Linear recalibration helpers
    // =========================================================================

    /// <summary>
    /// For each training mutation, clone the frozen graph, apply the mutation,
    /// and return the hop-1 MJ-direct DDG alongside the experimental value.
    /// No BFS convergence is run — this is a fast forward-only pass.
    /// </summary>
    private static List<(double mjDirect, double experimental)> GetTrainMjDirects(
        List<MutRecord> train,
        InMemoryGraph frozenGraph)
    {
        var results = new List<(double, double)>();
        foreach (var mut in train)
        {
            try
            {
                var graphCopy = frozenGraph.Clone();
                var agent = new MutationAgent(graphCopy);
                var (trace, _) = agent.ApplyMutationAndGetTrace(mut.ResidueId, mut.MutantType3);
                double mjDirect = trace.Where(t => t.hopDistance == 1).Sum(t => t.propagatedDelta);
                // Apply reference-state correction (same as RunOneMutationAsync)
                if (!string.IsNullOrEmpty(mut.WildType3))
                {
                    int nHop1 = trace.Count(t => t.hopDistance == 1);
                    mjDirect -= nHop1 * (MiyazawaJernigan.GetReferenceEnergy(mut.MutantType3)
                                       - MiyazawaJernigan.GetReferenceEnergy(mut.WildType3));
                }
                results.Add((mjDirect, mut.ExperimentalDdg));
            }
            catch
            {
                // Skip mutations that fail; they won't skew the fit
            }
        }
        return results;
    }

    /// <summary>
    /// Ordinary least-squares: fit experimental = slope * mjDirect + intercept
    /// on the provided (x, y) pairs.  Returns (1.0, 0.0) if n &lt; 2 or degenerate.
    /// </summary>
    private static (double slope, double intercept) FitOLS(
        List<(double x, double y)> pairs)
    {
        int n = pairs.Count;
        if (n < 2) return (1.0, 0.0);

        double mx = pairs.Average(p => p.x);
        double my = pairs.Average(p => p.y);
        double num = pairs.Sum(p => (p.x - mx) * (p.y - my));
        double den = pairs.Sum(p => (p.x - mx) * (p.x - mx));

        if (Math.Abs(den) < 1e-12) return (1.0, my - mx);
        double slope = num / den;
        double intercept = my - slope * mx;
        return (slope, intercept);
    }

    // =========================================================================
    // Correlation helpers
    // =========================================================================

    private static double PearsonR(IEnumerable<double> xs, IEnumerable<double> ys)
    {
        var xl = xs.ToList();
        var yl = ys.ToList();
        int n = Math.Min(xl.Count, yl.Count);
        if (n < 2) return double.NaN;

        double mx = xl.Take(n).Average();
        double my = yl.Take(n).Average();
        double num = 0, dxSq = 0, dySq = 0;
        for (int i = 0; i < n; i++)
        {
            double dx = xl[i] - mx, dy = yl[i] - my;
            num  += dx * dy;
            dxSq += dx * dx;
            dySq += dy * dy;
        }
        double denom = Math.Sqrt(dxSq * dySq);
        return denom < 1e-12 ? double.NaN : num / denom;
    }

    private static double SpearmanRho(IEnumerable<double> xs, IEnumerable<double> ys)
    {
        var xl = xs.ToList();
        var yl = ys.ToList();
        int n = Math.Min(xl.Count, yl.Count);
        if (n < 2) return double.NaN;
        return PearsonR(Ranks(xl, n), Ranks(yl, n));
    }

    private static IEnumerable<double> Ranks(List<double> vals, int n)
    {
        // Fractional (average) ranking
        var indexed = vals.Take(n).Select((v, i) => (v, i)).OrderBy(t => t.v).ToList();
        var ranks   = new double[n];
        int j = 0;
        while (j < n)
        {
            int k = j;
            while (k < n - 1 && indexed[k + 1].v == indexed[k].v) k++;
            double avgRank = (j + k + 2) / 2.0; // 1-based
            for (int m = j; m <= k; m++) ranks[indexed[m].i] = avgRank;
            j = k + 1;
        }
        return ranks;
    }

    // =========================================================================
    // Split generation / loading
    // =========================================================================

    private record MutRecord(
        string MutationId,
        string ResidueId,
        string WildType3,       // three-letter wildtype residue name from mutation code
        string MutantType3,
        double ExperimentalDdg,
        string Severity);    // "Low"|"Medium"|"High"

    private async Task<(List<MutRecord> train, List<MutRecord> val)> LoadOrGenerateSplitAsync(
        string familyName,
        string uniprotId,
        string csvPath,
        string splitFile)
    {
        // Try to load from cache first
        if (File.Exists(splitFile))
        {
            try
            {
                string json = await File.ReadAllTextAsync(splitFile);
                var cache = JsonSerializer.Deserialize<Dictionary<string, SplitCache>>(json) ?? new();
                if (cache.TryGetValue(familyName, out var cached))
                {
                    // Regenerate if cache predates WildType3 field (added for WT identity audit)
                    if (cached.Train.Count > 0 && string.IsNullOrEmpty(cached.Train[0].WildType3))
                    {
                        Console.WriteLine($" Split cache missing WildType3 — regenerating.");
                    }
                    else
                    {
                        Console.WriteLine($" Loaded split from {splitFile}.");
                        return (
                            cached.Train.Select(ToMutRecord).ToList(),
                            cached.Val.Select(ToMutRecord).ToList());
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($" Warning: could not read split cache ({ex.Message}), regenerating.");
            }
        }

        // Generate split from CSV
        var all = ReadMutationsFromCsv(csvPath, uniprotId);

        // Deduplicate by (MutationId) — same position+mutant from different experiments
        // should not inflate either train or val. Keep the first occurrence.
        var seen = new HashSet<string>();
        all = all.Where(m => seen.Add(m.MutationId)).ToList();

        Console.WriteLine($" Read {all.Count} mutations for {uniprotId} from CSV (after dedup).");

        if (all.Count == 0)
        {
            // Fallback: try reading from Neo4j
            Console.WriteLine(" Falling back to Neo4j family database...");
            all = await ReadMutationsFromNeo4jAsync(familyName, uniprotId);
            Console.WriteLine($" Read {all.Count} mutations from Neo4j.");
        }

        var (train, val) = StratifiedSplit(all, trainFraction: 0.80, seed: 42);
        await SaveSplitToCacheAsync(splitFile, familyName, train, val);
        Console.WriteLine($" Generated and saved split to {splitFile}.");
        return (train, val);
    }

    /// <summary>Stratified 80/20 split by severity (|ΔΔG| < 1 / 1–3 / > 3 kcal/mol).</summary>
    private static (List<MutRecord> train, List<MutRecord> val) StratifiedSplit(
        List<MutRecord> all, double trainFraction, int seed)
    {
        var rng   = new Random(seed);
        var train = new List<MutRecord>();
        var val   = new List<MutRecord>();

        foreach (var group in all.GroupBy(m => m.Severity))
        {
            var shuffled   = group.OrderBy(_ => rng.Next()).ToList();
            int trainCount = (int)Math.Round(shuffled.Count * trainFraction);
            train.AddRange(shuffled.Take(trainCount));
            val.AddRange(shuffled.Skip(trainCount));
        }
        return (train, val);
    }

    // =========================================================================
    // CSV reading
    // =========================================================================

    // AlphaFold/RCSB graphs always use chain "A" (single-chain structures).
    // The experimental CSV may use different chain labels (e.g., CI2 uses "I" in PDB).
    // This table maps each UniProt ID's CSV chain to the graph's actual chain ID.
    private static readonly Dictionary<string, string> GraphChain = new()
    {
        {"P00720", "A"},   // T4 lysozyme — 2LZM chain A
        {"P01053", "A"},   // CI2 — AlphaFold chain A (CSV uses "I")
        {"P00648", "A"},   // Barnase — AlphaFold chain A
    };

    // Offset to add to S2648 sequence positions to get the CIF/graph label_seq_id.
    // S2648 uses mature-protein (PDB) numbering; AlphaFold CIFs use UniProt canonical
    // numbering which includes the signal/propeptide prefix.
    // Empirically verified from the WT-audit mismatch patterns (Phase 9-D):
    //   CI2 P01053: mature chain starts at CIF position 2  → offset +1
    //   Barnase P00648: mature chain starts at CIF position 48 → offset +47
    //   T4 lysozyme P00720: RCSB 2LZM already uses PDB numbering → offset 0
    private static readonly Dictionary<string, int> SeqOffset = new()
    {
        {"P00720",  0},   // T4 lysozyme — RCSB, no offset
        {"P01053",  1},   // CI2 — 1-residue propeptide prefix in AlphaFold CIF
        {"P00648", 47},   // Barnase — 47-residue signal+propeptide in AlphaFold CIF
    };

    private static List<MutRecord> ReadMutationsFromCsv(string csvPath, string uniprotId)
    {
        if (!File.Exists(csvPath))
            return [];

        var records = new List<MutRecord>();
        using var reader = new StreamReader(csvPath);
        reader.ReadLine(); // skip header

        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var parts = line.Split(',');
            if (parts.Length < 7) continue;

            string uniprot  = parts[1].Trim();
            if (uniprot != uniprotId) continue;

            string chain        = parts[3].Trim();
            string mutationCode = parts[5].Trim();
            string ddgStr       = parts[6].Trim();

            if (!double.TryParse(ddgStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double ddg))
                continue;

            var m = Regex.Match(mutationCode, @"^([A-Z])(\d+)([A-Z])$", RegexOptions.IgnoreCase);
            if (!m.Success) continue;

            char wildChar  = char.ToUpperInvariant(m.Groups[1].Value[0]);
            int  seqPos    = int.Parse(m.Groups[2].Value)
                           + SeqOffset.GetValueOrDefault(uniprotId, 0);
            char mutantChar = char.ToUpperInvariant(m.Groups[3].Value[0]);
            if (!OneTo3.TryGetValue(wildChar,   out string? wildType3))   continue;
            if (!OneTo3.TryGetValue(mutantChar, out string? mutantType3)) continue;

            // Map CSV chain → graph chain. AlphaFold structures always use chain "A";
            // experimental databases (e.g., CI2 PDB "I") differ from the graph chain.
            string graphChain = GraphChain.TryGetValue(uniprotId, out string? gc) ? gc : chain;

            string residueId  = $"{uniprotId}_{graphChain}_{seqPos}";
            string mutationId = $"MUT_{uniprotId}_{graphChain}_{seqPos}_{mutantType3}";
            string severity   = Math.Abs(ddg) < 1.0 ? "Low"
                              : Math.Abs(ddg) <= 3.0 ? "Medium"
                              : "High";

            records.Add(new MutRecord(mutationId, residueId, wildType3, mutantType3, ddg, severity));
        }
        return records;
    }

    // =========================================================================
    // Neo4j fallback reading
    // =========================================================================

    private async Task<List<MutRecord>> ReadMutationsFromNeo4jAsync(string familyName, string uniprotId)
    {
        var records = new List<MutRecord>();
        try
        {
            await using var driver  = GraphDatabase.Driver(_uri, AuthTokens.Basic(_username, _password));
            await using var session = driver.AsyncSession(o => o.WithDatabase(familyName));
            await session.ExecuteReadAsync(async tx =>
            {
                var cursor = await tx.RunAsync(@"
                    MATCH (r:Residue)-[:HAS_MUTATION]->(m:Mutation)
                    WHERE m.uniprot = $u
                    RETURN m.id AS id, r.id AS residueId,
                           m.seqPos AS seqPos, m.mutantType AS mutantType,
                           m.ddg_kcal_mol AS ddg, m.mutatedChain AS chain",
                    new { u = uniprotId });

                while (await cursor.FetchAsync())
                {
                    double ddg = cursor.Current["ddg"].As<double>();
                    string severity = Math.Abs(ddg) < 1.0 ? "Low"
                                    : Math.Abs(ddg) <= 3.0 ? "Medium"
                                    : "High";
                    records.Add(new MutRecord(
                        cursor.Current["id"].As<string>(),
                        cursor.Current["residueId"].As<string>(),
                        "",   // WildType3 not stored in Neo4j fallback — WT audit skipped
                        cursor.Current["mutantType"].As<string>(),
                        ddg,
                        severity));
                }
            });
        }
        catch (Exception ex)
        {
            Console.WriteLine($" Neo4j fallback failed: {ex.Message}");
        }
        return records;
    }

    // =========================================================================
    // Split cache (JSON)
    // =========================================================================

    private sealed class SplitCacheEntry
    {
        public string MutationId    { get; set; } = "";
        public string ResidueId     { get; set; } = "";
        public string WildType3     { get; set; } = "";  // required for WT identity audit
        public string MutantType3   { get; set; } = "";
        public double ExperimentalDdg { get; set; }
        public string Severity      { get; set; } = "";
    }

    private sealed class SplitCache
    {
        public List<SplitCacheEntry> Train { get; set; } = [];
        public List<SplitCacheEntry> Val   { get; set; } = [];
    }

    private static MutRecord ToMutRecord(SplitCacheEntry e) =>
        new(e.MutationId, e.ResidueId, e.WildType3, e.MutantType3, e.ExperimentalDdg, e.Severity);

    private static async Task SaveSplitToCacheAsync(
        string splitFile, string familyName,
        List<MutRecord> train, List<MutRecord> val)
    {
        Dictionary<string, SplitCache> cache = new();
        if (File.Exists(splitFile))
        {
            try
            {
                cache = JsonSerializer.Deserialize<Dictionary<string, SplitCache>>(
                    await File.ReadAllTextAsync(splitFile)) ?? new();
            }
            catch { }
        }

        static SplitCacheEntry ToEntry(MutRecord r) => new()
        {
            MutationId     = r.MutationId,
            ResidueId      = r.ResidueId,
            WildType3      = r.WildType3,
            MutantType3    = r.MutantType3,
            ExperimentalDdg = r.ExperimentalDdg,
            Severity       = r.Severity
        };

        cache[familyName] = new SplitCache
        {
            Train = train.Select(ToEntry).ToList(),
            Val   = val.Select(ToEntry).ToList()
        };

        Directory.CreateDirectory(Path.GetDirectoryName(splitFile)!);
        await File.WriteAllTextAsync(splitFile,
            JsonSerializer.Serialize(cache, new JsonSerializerOptions { WriteIndented = true }));
    }
}

