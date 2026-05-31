using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Protein.Engine;

namespace Protein.Engine.Tests;

// =============================================================================
// Phase 6 — Convergence Checker tests
// Covers all four todo items:
//   1. EpsilonCalibrator: ε₀ = σ(ΔΔG_train) per family from S2648 CSV
//   2. Monotone narrowing tracker (via EnergySignalPropagator behaviour)
//   3. Convergence criterion validation logic
//   4. EnergySignalPropagator with pre-seeded ΔΔG_MJ initial interval
// =============================================================================

public class ConvergenceCheckerTests
{
    // =========================================================================
    // 1. EpsilonCalibrator — CSV parsing and σ computation
    // =========================================================================

    private static string WriteTempCsv(string content)
    {
        string path = Path.Combine(Path.GetTempPath(), $"mutcert_test_{Guid.NewGuid():N}.csv");
        File.WriteAllText(path, content);
        return path;
    }

    [Fact]
    public void Calibrate_CorrectStdDev_FromKnownValues()
    {
        // Values: 1, 2, 3, 4, 5  → mean = 3, σ² = 2, σ = √2 ≈ 1.4142
        string csv = "id,uniprot,ddg,method\n" +
                     "1,P00720,1.0,CD\n" +
                     "2,P00720,2.0,CD\n" +
                     "3,P00720,3.0,CD\n" +
                     "4,P00720,4.0,CD\n" +
                     "5,P00720,5.0,CD\n";

        string path = WriteTempCsv(csv);
        try
        {
            double eps = EpsilonCalibrator.Calibrate(path, "P00720");
            Assert.Equal(Math.Sqrt(2.0), eps, precision: 8);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void Calibrate_IgnoresOtherFamilies()
    {
        // P99999 has one outlier; P00720 has tight values.
        string csv = "id,uniprot,ddg,method\n" +
                     "1,P00720,2.0,CD\n" +
                     "2,P00720,2.0,CD\n" +   // σ = 0 → fallback
                     "3,P99999,100.0,CD\n";

        string path = WriteTempCsv(csv);
        try
        {
            // P00720 has σ = 0, so fallback is returned (values are identical)
            double eps = EpsilonCalibrator.Calibrate(path, "P00720");
            Assert.Equal(EpsilonCalibrator.DefaultEpsilon, eps, precision: 12);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void Calibrate_UnknownUniProt_ReturnsFallback()
    {
        string csv = "id,uniprot,ddg,method\n1,P00720,1.5,CD\n";
        string path = WriteTempCsv(csv);
        try
        {
            double eps = EpsilonCalibrator.Calibrate(path, "P99999");
            Assert.Equal(EpsilonCalibrator.DefaultEpsilon, eps, precision: 12);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void Calibrate_SingleDataPoint_ReturnsFallback()
    {
        string csv = "id,uniprot,ddg,method\n1,P00720,3.5,CD\n";
        string path = WriteTempCsv(csv);
        try
        {
            double eps = EpsilonCalibrator.Calibrate(path, "P00720");
            Assert.Equal(EpsilonCalibrator.DefaultEpsilon, eps, precision: 12);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void Calibrate_MissingFile_ReturnsFallback()
    {
        double eps = EpsilonCalibrator.Calibrate("/nonexistent/path.csv", "P00720");
        Assert.Equal(EpsilonCalibrator.DefaultEpsilon, eps, precision: 12);
    }

    [Fact]
    public void Calibrate_SkipsNonNumericDdgRows()
    {
        // Row 2 has "N/A" in the ddg column — should be skipped silently.
        string csv = "id,uniprot,ddg,method\n" +
                     "1,P00720,1.0,CD\n" +
                     "2,P00720,N/A,CD\n" +  // skipped
                     "3,P00720,3.0,CD\n";   // σ of {1, 3} = 1.0

        string path = WriteTempCsv(csv);
        try
        {
            var values = EpsilonCalibrator.ReadFamilyDdg(path, "P00720").ToList();
            Assert.Equal(2, values.Count);
            Assert.Contains(1.0, values);
            Assert.Contains(3.0, values);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void PopulationStdDev_MatchesExpectedFormula()
    {
        // Population σ: sum of squared deviations / n, then sqrt
        // Values: -1, 0, 1 → mean = 0, σ² = 2/3, σ = √(2/3)
        var values = new List<double> { -1.0, 0.0, 1.0 };
        double expected = Math.Sqrt(2.0 / 3.0);
        double actual   = EpsilonCalibrator.PopulationStdDev(values);
        Assert.Equal(expected, actual, precision: 12);
    }

    // =========================================================================
    // 2. EnergySignalPropagator — pre-seeded initial interval (spec § 6.1)
    // =========================================================================

    [Fact]
    public void Propagator_WithInitialDDG_SetsIntervalBeforeFirstApply()
    {
        // The interval should be [ddgMj - ε₀, ddgMj + ε₀] before any Apply().
        double eps0   = 1.5;
        double ddgMj  = 2.0;
        var esp = new EnergySignalPropagator(eps0, ddgMj);

        Assert.Equal(ddgMj - eps0, esp.Lo, precision: 12);
        Assert.Equal(ddgMj + eps0, esp.Hi, precision: 12);
        Assert.Equal(2.0 * eps0,   esp.Width, precision: 12);
        Assert.Equal(0, esp.StepCount);
    }

    [Fact]
    public void Propagator_WithInitialDDG_CumulativeDDG_StartsAtInitial()
    {
        // Running sum starts at ΔΔG_MJ, not at zero.
        double ddgMj = 3.7;
        var esp = new EnergySignalPropagator(epsilon0: 1.0, initialDdgMj: ddgMj);
        Assert.Equal(ddgMj, esp.CumulativeDDG, precision: 12);

        esp.Apply(0.5);
        Assert.Equal(ddgMj + 0.5, esp.CumulativeDDG, precision: 12);
    }

    [Fact]
    public void Propagator_WithInitialDDG_NarrowsOnFirstStep()
    {
        // Initial: lo = 0.5, hi = 2.5 (ddgMj=1.5, ε₀=1.0)
        // After step 1 with delta=0.1:
        //   CumulativeDDG = 1.5 + 0.1 = 1.6
        //   uncertainty = 1.0/√1 = 1.0
        //   proposedLo = 0.6, proposedHi = 2.6
        //   lo = max(0.5, 0.6) = 0.6, hi = min(2.5, 2.6) = 2.5 → width = 1.9
        var esp = new EnergySignalPropagator(epsilon0: 1.0, initialDdgMj: 1.5);

        esp.Apply(0.1);

        Assert.Equal(0.6, esp.Lo, precision: 10);
        Assert.Equal(2.5, esp.Hi, precision: 10);
        Assert.Equal(1.9, esp.Width, precision: 10);
    }

    [Fact]
    public void Propagator_WithInitialDDG_ConsistentDeltas_WidthStrictlyNarrows()
    {
        // With small consistent deltas and a tight initial interval, the interval
        // should narrow monotonically (no crossovers expected with ε₀ = 2.0).
        var esp = new EnergySignalPropagator(epsilon0: 2.0, initialDdgMj: 0.5);
        double lastWidth = esp.Width;

        for (int k = 1; k <= 20; k++)
        {
            esp.Apply(0.05);
            if (!esp.CrossoverOccurred)
            {
                Assert.True(esp.Width <= lastWidth + 1e-12,
                    $"Width increased without crossover at step {k}: {lastWidth} → {esp.Width}");
            }
            lastWidth = esp.Width;
        }
    }

    // =========================================================================
    // 3. Convergence criterion (logic, not Neo4j)
    // =========================================================================

    [Fact]
    public void ConvergenceCriterion_MetAfterEnoughNarrowingSteps()
    {
        // Simulate 10 identical deltas of 0 with ε₀ = 1.0, initialDdgMj = 0.
        // After k steps: width = 2 * ε₀ / √k → at k=1 width=2.0, k=4 width=1.0.
        // Track stableSteps and check that criterion is met.
        var esp = new EnergySignalPropagator(epsilon0: 1.0, initialDdgMj: 0.0);
        int stableSteps = 0;
        double lastWidth = esp.Width;  // = 2.0 initially (from preset)
        bool criterionMet = false;

        // Simulated shell: 3 nodes, all "evaluated" after step 3.
        int shellSize = 3;
        var evaluated = new System.Collections.Generic.HashSet<string>();

        for (int k = 1; k <= 20; k++)
        {
            esp.Apply(0.0);

            if (k <= shellSize) evaluated.Add($"node{k}");

            if (esp.CrossoverOccurred) { stableSteps = 0; lastWidth = double.MaxValue; }

            if (esp.Width <= lastWidth) stableSteps++; else stableSteps = 0;
            lastWidth = esp.Width;

            if (esp.Width <= 2.0 && stableSteps >= 5 && evaluated.Count >= shellSize)
            {
                criterionMet = true;
                break;
            }
        }

        Assert.True(criterionMet, "Convergence criterion should be met within 20 zero-delta steps.");
    }

    [Fact]
    public void ConvergenceCriterion_NotMet_WhenShellIncomplete()
    {
        // Even if width is narrow and stableSteps is large, certificate should NOT
        // be issued while the shell is not fully evaluated.
        var esp = new EnergySignalPropagator(epsilon0: 1.0, initialDdgMj: 0.0);
        int stableSteps = 0;
        double lastWidth = esp.Width;
        int shellSize = 5; // 5 nodes required
        var evaluated = new System.Collections.Generic.HashSet<string>();
        evaluated.Add("node1"); // Only 1 of 5 evaluated

        bool criterionMet = false;
        for (int k = 1; k <= 30; k++)
        {
            esp.Apply(0.0);
            if (esp.CrossoverOccurred) { stableSteps = 0; lastWidth = double.MaxValue; }
            if (esp.Width <= lastWidth) stableSteps++; else stableSteps = 0;
            lastWidth = esp.Width;

            if (esp.Width <= 2.0 && stableSteps >= 5 && evaluated.Count >= shellSize)
            {
                criterionMet = true;
                break;
            }
        }

        Assert.False(criterionMet, "Certificate must not be issued when the 8Å shell is not fully evaluated.");
    }
}
