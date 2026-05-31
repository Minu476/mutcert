using System;
using System.Collections.Generic;

namespace Protein.Engine;

/// <summary>
/// Provides high-performance, O(1) lookups for the Miyazawa-Jernigan (MJ) 1996 contact potentials
/// (AAindex entry MIYS960101). Energies are converted from statistical units (k_B T) to physical kcal/mol at 298K.
/// </summary>
public static class MiyazawaJernigan
{
    private const double ConversionFactor = 0.5922; // Convert k_B T to kcal/mol at 298 K
    private const string AminoAcidOrder = "ARNDCQEGHILKMFPSTWYV";

    private static readonly Dictionary<string, int> ThreeToOneMap = new(StringComparer.OrdinalIgnoreCase)
    {
        { "ALA", 0 }, { "ARG", 1 }, { "ASN", 2 }, { "ASP", 3 }, { "CYS", 4 },
        { "GLN", 5 }, { "GLU", 6 }, { "GLY", 7 }, { "HIS", 8 }, { "ILE", 9 },
        { "LEU", 10 }, { "LYS", 11 }, { "MET", 12 }, { "PHE", 13 }, { "PRO", 14 },
        { "SER", 15 }, { "THR", 16 }, { "TRP", 17 }, { "TYR", 18 }, { "VAL", 19 }
    };

    private static readonly double[,] EnergyMatrix;

    /// <summary>
    /// Swiss-Prot background amino acid frequencies, same order as AminoAcidOrder ("ARNDCQEGHILKMFPSTWYV").
    /// Used to compute the mean-field unfolded-state reference energy for each residue type.
    /// Source: UniProt Swiss-Prot release 2024 amino acid composition statistics.
    /// </summary>
    private static readonly double[] BackgroundFrequency =
    [
        0.0777, // A
        0.0530, // R
        0.0430, // N
        0.0535, // D
        0.0184, // C
        0.0406, // Q
        0.0635, // E
        0.0688, // G
        0.0228, // H
        0.0565, // I
        0.0965, // L
        0.0581, // K
        0.0241, // M
        0.0390, // F
        0.0471, // P
        0.0684, // S
        0.0568, // T
        0.0112, // W
        0.0290, // Y
        0.0664  // V
    ];

    /// <summary>
    /// Pre-computed mean-field reference energy for each residue type:
    /// ref(X) = Σ_j f(j) × E(X, j) where f(j) is the Swiss-Prot background frequency.
    /// Used for unfolded-state reference-state subtraction.
    /// </summary>
    private static readonly double[] ReferenceEnergy;

    static MiyazawaJernigan()
    {
        // 20x20 Symmetric contact potential matrix in k_B T units (MIYS960101 from AAindex database)
        double[][] lowerTriangle = [
            [-2.72], // A
            [-1.83, -1.55], // R
            [-1.84, -1.64, -1.68], // N
            [-1.70, -2.29, -1.68, -1.21], // D
            [-3.57, -2.57, -2.59, -2.41, -5.44], // C
            [-1.89, -1.80, -1.71, -1.46, -2.85, -1.54], // Q
            [-1.51, -2.27, -1.51, -1.02, -2.27, -1.42, -0.91], // E
            [-2.31, -1.72, -1.74, -1.59, -3.16, -1.66, -1.22, -2.24], // G
            [-2.41, -2.16, -2.08, -2.32, -3.60, -1.98, -2.15, -2.15, -3.05], // H
            [-4.58, -3.63, -3.24, -3.17, -5.50, -3.67, -3.27, -3.78, -4.14, -6.54], // I
            [-4.91, -4.03, -3.74, -3.40, -5.83, -4.04, -3.59, -4.16, -4.54, -7.04, -7.37], // L
            [-1.31, -0.59, -1.21, -1.68, -1.95, -1.29, -1.80, -1.15, -1.35, -3.01, -3.37, -0.12], // K
            [-3.94, -3.12, -2.95, -2.57, -4.99, -3.30, -2.89, -3.39, -3.98, -6.02, -6.41, -2.48, -5.46], // M
            [-4.81, -3.98, -3.75, -3.48, -5.80, -4.10, -3.56, -4.13, -4.77, -6.84, -7.28, -3.36, -6.56, -7.26], // F
            [-2.03, -1.70, -1.53, -1.33, -3.07, -1.73, -1.26, -1.87, -2.25, -3.76, -4.20, -0.97, -3.45, -4.25, -1.75], // P
            [-2.01, -1.62, -1.58, -1.63, -2.86, -1.49, -1.48, -1.82, -2.11, -3.52, -3.92, -1.05, -3.03, -4.02, -1.57, -1.67], // S
            [-2.32, -1.90, -1.88, -1.80, -3.11, -1.90, -1.74, -2.08, -2.42, -4.03, -4.34, -1.31, -3.51, -4.28, -1.90, -1.96, -2.12], // T
            [-3.82, -3.41, -3.07, -2.84, -4.95, -3.11, -2.99, -3.42, -3.98, -5.78, -6.14, -2.69, -5.55, -6.16, -3.73, -2.99, -3.22, -5.06], // W
            [-3.36, -3.16, -2.76, -2.76, -4.16, -2.97, -2.79, -3.01, -3.52, -5.25, -5.67, -2.60, -4.91, -5.66, -3.19, -2.78, -3.01, -4.66, -4.17], // Y
            [-4.04, -3.07, -2.83, -2.48, -4.96, -3.07, -2.67, -3.38, -3.58, -6.05, -6.48, -2.49, -5.32, -6.29, -3.32, -3.05, -3.46, -5.18, -4.62, -5.52]  // V
        ];

        EnergyMatrix = new double[20, 20];
        for (int i = 0; i < 20; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double rawEnergy = lowerTriangle[i][j];
                double energyKcal = rawEnergy * ConversionFactor;
                EnergyMatrix[i, j] = energyKcal;
                EnergyMatrix[j, i] = energyKcal;
            }
        }

        // Pre-compute reference energies: ref(X) = Σ_j f(j) × E(X, j)
        ReferenceEnergy = new double[20];
        for (int i = 0; i < 20; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < 20; j++)
                sum += BackgroundFrequency[j] * EnergyMatrix[i, j];
            ReferenceEnergy[i] = sum;
        }
    }

    /// <summary>
    /// Returns the mean-field reference energy for a residue type in the unfolded (random-coil) state.
    /// ref(X) = Σ_j f(j) × E(X, j), where f(j) is the Swiss-Prot background frequency.
    /// Used to subtract the unfolded-state contribution from ΔΔG_MJ, bringing its scale
    /// in line with experimental ΔΔG values.
    /// </summary>
    public static double GetReferenceEnergy(string residueType)
    {
        int idx = ResolveIndex(residueType);
        return ReferenceEnergy[idx];
    }

    /// <summary>
    /// Checks if a 3-letter or 1-letter residue represents a hydrophobic amino acid.
    /// Hydrophobic amino acids are ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, TYR.
    /// </summary>
    public static bool IsHydrophobic(string residueCode)
    {
        string formatted = residueCode.Trim().ToUpperInvariant();
        if (formatted.Length == 3)
        {
            return formatted is "ALA" or "VAL" or "LEU" or "ILE" or "MET" or "PHE" or "TRP" or "PRO" or "TYR";
        }
        else if (formatted.Length == 1)
        {
            return formatted is "A" or "V" or "L" or "I" or "M" or "F" or "W" or "P" or "Y";
        }
        return false;
    }

    /// <summary>
    /// Gets the charge of a standard residue (+1.0, -1.0, +0.5, or 0.0).
    /// </summary>
    public static double GetCharge(string residueCode)
    {
        string formatted = residueCode.Trim().ToUpperInvariant();
        return formatted switch
        {
            "ASP" or "D" => -1.0,
            "GLU" or "E" => -1.0,
            "LYS" or "K" => 1.0,
            "ARG" or "R" => 1.0,
            "HIS" or "H" => 0.5,
            _ => 0.0
        };
    }

    /// <summary>
    /// Retrieves the contact energy between two amino acids in kcal/mol.
    /// </summary>
    /// <param name="aa1">Three-letter or single-letter code for residue 1 (e.g. MET, M)</param>
    /// <param name="aa2">Three-letter or single-letter code for residue 2 (e.g. CYS, C)</param>
    public static double GetEnergy(string aa1, string aa2)
    {
        int idx1 = ResolveIndex(aa1);
        int idx2 = ResolveIndex(aa2);

        if (idx1 < 0 || idx2 < 0)
        {
            throw new ArgumentException($"Unsupported or invalid amino acid combination: '{aa1}' and '{aa2}'");
        }

        return EnergyMatrix[idx1, idx2];
    }

    private static int ResolveIndex(string aa)
    {
        string formatted = aa.Trim().ToUpperInvariant();
        if (formatted == "MSE") formatted = "MET";
        if (formatted == "ASX") formatted = "ASP";
        if (formatted == "GLX") formatted = "GLU";
        if (formatted == "SEC") formatted = "CYS";
        if (formatted.Length == 3)
        {
            return ThreeToOneMap.TryGetValue(formatted, out int idx) ? idx : ThreeToOneMap["ALA"]; // Fallback to ALA
        }
        else if (formatted.Length == 1)
        {
            int idx = AminoAcidOrder.IndexOf(formatted[0]);
            return idx >= 0 ? idx : 0; // Fallback to index 0 (ALA)
        }
        return 0; // Fallback to index 0 (ALA)
    }
}
