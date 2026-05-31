using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Protein.Engine;

/// <summary>
/// Computes ε₀ = σ(ΔΔG_train) per protein family from the S2648 / ThermoMutDB CSV.
///
/// Per the MutCert spec (v3.3):
///   ε₀ = σ(ΔΔG_{train, family})
///
/// This calibration parameter sets the initial interval half-width used by
/// <see cref="EnergySignalPropagator"/> before BFS propagation begins.
///
/// Usage:
///   double eps0 = EpsilonCalibrator.Calibrate("data/s2648/s2648.csv", "P00720");
/// </summary>
public static class EpsilonCalibrator
{
    /// <summary>
    /// Default ε₀ returned when fewer than 2 data points are available for a family.
    /// Matches the Phase 1 hardcoded value used in earlier milestones.
    /// </summary>
    public const double DefaultEpsilon = 1.5;

    /// <summary>
    /// Reads <paramref name="csvPath"/>, filters rows whose <c>uniprot</c> column matches
    /// <paramref name="uniprotId"/>, and returns the population standard deviation of the
    /// <c>ddg</c> column.
    ///
    /// Returns <paramref name="fallback"/> (default 1.5) when:
    /// - the file does not exist
    /// - fewer than 2 data points are found for the family
    /// - the standard deviation is zero (degenerate data)
    /// </summary>
    public static double Calibrate(
        string csvPath,
        string uniprotId,
        double fallback = DefaultEpsilon)
    {
        var ddgValues = ReadFamilyDdg(csvPath, uniprotId).ToList();

        if (ddgValues.Count < 2)
        {
            Console.WriteLine(
                $"[EpsilonCalibrator] Insufficient data for {uniprotId} " +
                $"({ddgValues.Count} record(s)). Using fallback ε₀ = {fallback}.");
            return fallback;
        }

        double sigma = PopulationStdDev(ddgValues);

        if (sigma == 0.0)
        {
            Console.WriteLine(
                $"[EpsilonCalibrator] All ΔΔG values identical for {uniprotId}. " +
                $"Using fallback ε₀ = {fallback}.");
            return fallback;
        }

        Console.WriteLine(
            $"[EpsilonCalibrator] {uniprotId}: n = {ddgValues.Count}, " +
            $"σ(ΔΔG) = {sigma:F4} kcal/mol → ε₀ = {sigma:F4}");

        return sigma;
    }

    // -------------------------------------------------------------------------
    // Helpers (public for testability)
    // -------------------------------------------------------------------------

    /// <summary>
    /// Streams <c>ddg</c> values from the CSV for the given UniProt ID.
    /// Skips rows with non-numeric or missing ddg values.
    /// </summary>
    public static IEnumerable<double> ReadFamilyDdg(string csvPath, string uniprotId)
    {
        if (!File.Exists(csvPath))
            yield break;

        bool firstLine = true;
        int uniprotCol = -1;
        int ddgCol     = -1;

        foreach (string line in File.ReadLines(csvPath))
        {
            // Skip blank lines
            if (string.IsNullOrWhiteSpace(line))
                continue;

            string[] cols = line.Split(',');

            if (firstLine)
            {
                for (int i = 0; i < cols.Length; i++)
                {
                    string header = cols[i].Trim().ToLowerInvariant();
                    if (header == "uniprot") uniprotCol = i;
                    if (header == "ddg")     ddgCol     = i;
                }
                firstLine = false;
                continue;
            }

            if (uniprotCol < 0 || ddgCol < 0)
                continue;

            if (cols.Length <= Math.Max(uniprotCol, ddgCol))
                continue;

            if (!string.Equals(cols[uniprotCol].Trim(), uniprotId, StringComparison.OrdinalIgnoreCase))
                continue;

            string rawDdg = cols[ddgCol].Trim();
            if (double.TryParse(rawDdg, NumberStyles.Float, CultureInfo.InvariantCulture, out double ddg))
                yield return ddg;
        }
    }

    /// <summary>
    /// Population standard deviation: √( Σ(xᵢ − μ)² / n ).
    /// Returns 0.0 for lists with fewer than 2 elements (already checked by caller).
    /// </summary>
    public static double PopulationStdDev(IList<double> values)
    {
        if (values.Count < 2) return 0.0;
        double mean     = values.Average();
        double variance = values.Sum(v => (v - mean) * (v - mean)) / values.Count;
        return Math.Sqrt(variance);
    }

    /// <summary>
    /// Calibrates ε₀ directly from an in-memory sequence of ΔΔG values
    /// (typically the training split). No test-set data should be passed.
    /// </summary>
    public static double CalibrateFromValues(
        IEnumerable<double> trainDdg,
        string uniprotId,
        double fallback = DefaultEpsilon)
    {
        var values = trainDdg.ToList();
        if (values.Count < 2)
        {
            Console.WriteLine(
                $"[EpsilonCalibrator] Insufficient training data for {uniprotId} " +
                $"({values.Count} records). Using fallback ε₀ = {fallback}.");
            return fallback;
        }

        double sigma = PopulationStdDev(values);
        if (sigma == 0.0)
        {
            Console.WriteLine(
                $"[EpsilonCalibrator] All training ΔΔG values identical for {uniprotId}. " +
                $"Using fallback ε₀ = {fallback}.");
            return fallback;
        }

        Console.WriteLine(
            $"[EpsilonCalibrator] {uniprotId}: n_train = {values.Count}, " +
            $"σ(ΔΔG_train) = {sigma:F4} kcal/mol → ε₀ = {sigma:F4}");
        return sigma;
    }
}
