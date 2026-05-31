using System.Collections.Generic;
using System.Numerics;
using Protein.Engine;

namespace Protein.Engine.Tests;

/// <summary>
/// Unit tests for the three BrainCore guard conditions.
/// Each test builds a minimal synthetic residue list that represents
/// a known violation, verifies the guard rejects it (returns false),
/// and also verifies that a structurally valid list is accepted (returns true).
/// </summary>
public class BrainCoreTests
{
    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------

    /// <summary>
    /// Creates a minimal backbone residue at a given CA position.
    /// N and C are placed 1.0 Å either side of CA along the X-axis
    /// so backbone planarity can also be tested independently.
    /// </summary>
    private static Residue MakeResidue(
        string name,
        int seqPos,
        Vector3 ca,
        string chainId = "A",
        float charge = 0f)
    {
        var r = new Residue
        {
            ResidueName = name,
            SeqPos      = seqPos,
            ChainId     = chainId,
            PLDDT       = 90.0,
            CA          = ca,
            N           = ca - new Vector3(1.0f, 0f, 0f),
            C           = ca + new Vector3(1.0f, 0f, 0f),
            O           = ca + new Vector3(1.0f, 1.0f, 0f),
        };
        r.ComputeSideChainCentroid();
        return r;
    }

    // -------------------------------------------------------------------------
    // 1. Steric Clash Detection
    // -------------------------------------------------------------------------

    [Fact]
    public void StericClash_Violation_ReturnsFalse()
    {
        // Two non-adjacent residues (seqPos 1 and 10) whose CA atoms are
        // placed at the same coordinate → distance = 0, well below the 3.0 Å threshold.
        var clash = new List<Residue>
        {
            MakeResidue("ALA", 1,  new Vector3(0f, 0f, 0f)),
            MakeResidue("ALA", 10, new Vector3(0f, 0f, 0f))  // exact clash
        };

        bool result = BrainCore.ValidateStericClash(clash);

        Assert.False(result, "Expected steric clash to be detected (CA-CA distance = 0 Å).");
    }

    [Fact]
    public void StericClash_BorderlineViolation_ReturnsFalse()
    {
        // CA-CA = 2.5 Å — below the 3.0 Å minimum for non-adjacent residues.
        var residues = new List<Residue>
        {
            MakeResidue("ALA", 1,  new Vector3(0f, 0f, 0f)),
            MakeResidue("GLY", 5,  new Vector3(2.5f, 0f, 0f))
        };

        Assert.False(BrainCore.ValidateStericClash(residues),
            "CA-CA distance of 2.5 Å (< 3.0 Å) must be flagged as a steric clash.");
    }

    [Fact]
    public void StericClash_AdjacentResidues_NotFlagged()
    {
        // Adjacent residues (seqPos 1 and 2) are exempted even at CA-CA = 0.
        var residues = new List<Residue>
        {
            MakeResidue("ALA", 1, new Vector3(0f, 0f, 0f)),
            MakeResidue("ALA", 2, new Vector3(0f, 0f, 0f))
        };

        Assert.True(BrainCore.ValidateStericClash(residues),
            "Adjacent residues (|seqPos| ≤ 1) must be exempt from steric clash check.");
    }

    [Fact]
    public void StericClash_ValidStructure_ReturnsTrue()
    {
        // Three consecutive residues spaced 3.8 Å apart — a normal backbone Cα-Cα distance.
        var residues = new List<Residue>
        {
            MakeResidue("ALA", 1, new Vector3(0f,   0f, 0f)),
            MakeResidue("VAL", 2, new Vector3(3.8f, 0f, 0f)),
            MakeResidue("LEU", 3, new Vector3(7.6f, 0f, 0f))
        };

        Assert.True(BrainCore.ValidateStericClash(residues),
            "Residues spaced ≥ 3.0 Å apart must pass the steric clash check.");
    }

    // -------------------------------------------------------------------------
    // 2. Bond Planarity (backbone C–N peptide bond length)
    // -------------------------------------------------------------------------

    /// <summary>
    /// Builds a pair of consecutive residues where the C–N bond length is
    /// explicitly controlled by placing the atoms.
    /// </summary>
    private static List<Residue> MakePeptidePair(float cToNDistance)
    {
        // Residue 1: C at origin
        var r1 = new Residue
        {
            ResidueName = "ALA", SeqPos = 1, ChainId = "A",
            CA = new Vector3(-1f, 0f, 0f),
            N  = new Vector3(-2f, 0f, 0f),
            C  = new Vector3(0f,  0f, 0f),   // C at origin
            O  = new Vector3(0f,  1f, 0f),
        };
        r1.ComputeSideChainCentroid();

        // Residue 2: N placed exactly cToNDistance away from r1.C along X-axis
        var r2 = new Residue
        {
            ResidueName = "GLY", SeqPos = 2, ChainId = "A",
            N  = new Vector3(cToNDistance, 0f, 0f),
            CA = new Vector3(cToNDistance + 1f, 0f, 0f),
            C  = new Vector3(cToNDistance + 2f, 0f, 0f),
            O  = new Vector3(cToNDistance + 2f, 1f, 0f),
        };
        r2.ComputeSideChainCentroid();

        return [r1, r2];
    }

    [Fact]
    public void BackbonePlanarity_CNDistanceTooShort_ReturnsFalse()
    {
        // C–N = 0.9 Å, below the 1.17 Å minimum
        var residues = MakePeptidePair(cToNDistance: 0.9f);
        Assert.False(BrainCore.ValidateBackbonePlanarity(residues),
            "C–N distance of 0.9 Å (< 1.17 Å) must fail backbone planarity check.");
    }

    [Fact]
    public void BackbonePlanarity_CNDistanceTooLong_ReturnsFalse()
    {
        // C–N = 1.8 Å, above the 1.47 Å maximum
        var residues = MakePeptidePair(cToNDistance: 1.8f);
        Assert.False(BrainCore.ValidateBackbonePlanarity(residues),
            "C–N distance of 1.8 Å (> 1.47 Å) must fail backbone planarity check.");
    }

    [Fact]
    public void BackbonePlanarity_NormalCNDistance_ReturnsTrue()
    {
        // C–N = 1.32 Å — the ideal peptide bond length
        var residues = MakePeptidePair(cToNDistance: 1.32f);
        Assert.True(BrainCore.ValidateBackbonePlanarity(residues),
            "C–N distance of 1.32 Å (ideal peptide bond) must pass backbone planarity check.");
    }

    [Fact]
    public void BackbonePlanarity_BoundaryValues_PassesAtEdges()
    {
        // Exactly at the minimum (1.17) and maximum (1.47) — should pass
        var atMin = MakePeptidePair(cToNDistance: 1.17f);
        var atMax = MakePeptidePair(cToNDistance: 1.47f);

        Assert.True(BrainCore.ValidateBackbonePlanarity(atMin),
            "C–N = 1.17 Å (boundary min) must pass.");
        Assert.True(BrainCore.ValidateBackbonePlanarity(atMax),
            "C–N = 1.47 Å (boundary max) must pass.");
    }

    // -------------------------------------------------------------------------
    // 3. Local Charge Neutrality
    // -------------------------------------------------------------------------

    /// <summary>
    /// Builds a list of charged residues all within 10 Å of the origin.
    /// </summary>
    private static List<Residue> MakeChargedResidues(
        string[] residueNames,
        float spreadRadius = 2.0f)
    {
        var residues = new List<Residue>();
        for (int i = 0; i < residueNames.Length; i++)
        {
            // Spread residues in a small sphere within the test radius
            float x = spreadRadius * (i % 3) / 3f;
            float y = spreadRadius * (i / 3) / 3f;
            var r = MakeResidue(residueNames[i], i + 1, new Vector3(x, y, 0f));
            residues.Add(r);
        }
        return residues;
    }

    [Fact]
    public void LocalChargeNeutrality_ExcessivePositiveCharge_ReturnsFalse()
    {
        // 6× Arg (+1 each) = total charge +6.0, exceeds the ±5.0 threshold → FAIL
        var residues = MakeChargedResidues(["ARG", "ARG", "ARG", "ARG", "ARG", "ARG"]);
        var center = Vector3.Zero;

        Assert.False(BrainCore.ValidateLocalChargeNeutrality(residues, center, radius: 10.0f),
            "Total charge of +6.0 must fail the local charge neutrality check (|q| > 5.0).");
    }

    [Fact]
    public void LocalChargeNeutrality_ExcessiveNegativeCharge_ReturnsFalse()
    {
        // 6× Asp (-1 each) = total charge -6.0 → FAIL
        var residues = MakeChargedResidues(["ASP", "ASP", "ASP", "ASP", "ASP", "ASP"]);

        Assert.False(BrainCore.ValidateLocalChargeNeutrality(residues, Vector3.Zero, radius: 10.0f),
            "Total charge of -6.0 must fail the local charge neutrality check (|q| > 5.0).");
    }

    [Fact]
    public void LocalChargeNeutrality_NeutralMix_ReturnsTrue()
    {
        // 3× Arg (+3) + 3× Asp (-3) = total charge 0.0 → PASS
        var residues = MakeChargedResidues(["ARG", "ARG", "ARG", "ASP", "ASP", "ASP"]);

        Assert.True(BrainCore.ValidateLocalChargeNeutrality(residues, Vector3.Zero, radius: 10.0f),
            "Balanced +3 / -3 charge must pass the local charge neutrality check.");
    }

    [Fact]
    public void LocalChargeNeutrality_AllNeutralResidues_ReturnsTrue()
    {
        // Only Ala (charge 0) — total = 0 → PASS
        var residues = MakeChargedResidues(["ALA", "ALA", "ALA", "GLY", "VAL", "LEU"]);

        Assert.True(BrainCore.ValidateLocalChargeNeutrality(residues, Vector3.Zero, radius: 10.0f),
            "All-neutral residues must pass the local charge neutrality check.");
    }

    [Fact]
    public void LocalChargeNeutrality_ChargedResidueBeyondRadius_Ignored()
    {
        // 6× Arg placed 50 Å from center — outside the 10 Å radius → ignored → PASS
        var residues = new List<Residue>();
        for (int i = 0; i < 6; i++)
            residues.Add(MakeResidue("ARG", i + 1, new Vector3(50f + i, 0f, 0f)));

        Assert.True(BrainCore.ValidateLocalChargeNeutrality(residues, Vector3.Zero, radius: 10.0f),
            "Charged residues beyond the search radius must be ignored.");
    }
}
