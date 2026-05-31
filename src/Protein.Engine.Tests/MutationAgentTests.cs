using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Protein.Engine;

namespace Protein.Engine.Tests;

// =============================================================================
// Phase 5 — Mutation Agent tests
// Covers all four todo items:
//   1. Residue node swap logic (wildtype → mutant)
//   2. Edge re-evaluation per rules table (create / delete edges on swap)
//   3. BFS energy propagation within 8 Å shell (trace is non-empty, hop-1 has
//      the correct energy delta when the mutation changes MJ energies)
//   4. EnergySignalPropagator: update [lo, hi] interval after each BFS step
// =============================================================================

public class MutationAgentTests
{
    // -------------------------------------------------------------------------
    // Graph-building helpers
    // -------------------------------------------------------------------------

    private static Residue MakeRes(
        string name, int seqPos, Vector3 ca,
        Vector3? cb = null, string chain = "A")
    {
        var r = new Residue
        {
            ResidueName = name,
            SeqPos      = seqPos,
            ChainId     = chain,
            PLDDT       = 90.0,
            CA = ca,
            N  = ca - new Vector3(1.0f, 0f, 0f),
            C  = ca + new Vector3(1.0f, 0f, 0f),
            O  = ca + new Vector3(1.0f, 1.0f, 0f),
            CB = cb
        };
        r.ComputeSideChainCentroid();
        return r;
    }

    /// <summary>
    /// Builds a minimal 3-residue graph (chain A, seqPos 1-2-3) with a peptide
    /// backbone and an H-bond from residue 2 (the mutation target) to residue 3,
    /// where residue 2 is the donor.
    /// </summary>
    private static (InMemoryGraph graph, string id1, string id2, string id3)
        BuildProGraph()
    {
        const string uid = "TEST";
        string id1 = $"{uid}_A_1";
        string id2 = $"{uid}_A_2";
        string id3 = $"{uid}_A_3";

        var graph = new InMemoryGraph();
        graph.Nodes[id1] = MakeRes("ALA", 1, new Vector3(0f,   0f, 0f));
        graph.Nodes[id2] = MakeRes("LYS", 2, new Vector3(3.8f, 0f, 0f));
        graph.Nodes[id3] = MakeRes("ASP", 3, new Vector3(7.6f, 0f, 0f));

        // Peptide bonds (backbone)
        graph.Edges.Add(new Edge { FromId = id1, ToId = id2, Type = EdgeType.Peptide, EnergyKcal = 0 });
        graph.Edges.Add(new Edge { FromId = id2, ToId = id3, Type = EdgeType.Peptide, EnergyKcal = 0 });

        // H-bond: residue 2 (LYS) is donor, residue 3 (ASP) is acceptor
        graph.Edges.Add(new Edge
        {
            FromId     = id2,
            ToId       = id3,
            Type       = EdgeType.HBond,
            EnergyKcal = -0.5,
            Donor      = id2,
            Acceptor   = id3
        });

        graph.BuildAdjacencyList();
        return (graph, id1, id2, id3);
    }

    /// <summary>
    /// Builds a 2-residue graph where both residues are hydrophobic,
    /// their CB atoms are 4 Å apart, and they are connected by a
    /// HYDROPHOBIC_CONTACT edge plus a PEPTIDE bond.
    /// Residue 1 (GLU, non-hydrophobic) is the mutation target;
    /// residue 2 (ILE) is the neighbour.
    /// </summary>
    private static (InMemoryGraph graph, string targetId, string partnerId)
        BuildHydrophobicGraph()
    {
        const string uid = "TEST";
        string targetId  = $"{uid}_A_1";
        string partnerId = $"{uid}_A_5"; // non-adjacent seqPos so no adjacency exemption

        var targetRes  = MakeRes("GLU", 1, new Vector3(0f,   0f, 0f), cb: new Vector3(0.5f, 0f, 0f));
        var partnerRes = MakeRes("ILE", 5, new Vector3(4.0f, 0f, 0f), cb: new Vector3(4.5f, 0f, 0f));

        var graph = new InMemoryGraph();
        graph.Nodes[targetId]  = targetRes;
        graph.Nodes[partnerId] = partnerRes;

        // No hydrophobic edge initially (GLU is not hydrophobic)
        // Only electrostatic edge (GLU has charge -1, ILE has 0 → no electrostatic)
        // Just a peptide for connectivity
        graph.Edges.Add(new Edge { FromId = targetId, ToId = partnerId, Type = EdgeType.Peptide, EnergyKcal = 0 });

        graph.BuildAdjacencyList();
        return (graph, targetId, partnerId);
    }

    /// <summary>
    /// Builds a 2-residue graph: charged LYS target and charged ASP partner
    /// connected by an ELECTROSTATIC edge, separated by 5 Å at the side-chain centroids.
    /// </summary>
    private static (InMemoryGraph graph, string targetId, string partnerId)
        BuildElectrostaticGraph()
    {
        const string uid = "TEST";
        string targetId  = $"{uid}_A_1";
        string partnerId = $"{uid}_A_5";

        // Side-chain centroids default to CA when no side-chain atoms exist
        var targetRes  = MakeRes("LYS", 1, new Vector3(0f,   0f, 0f));
        var partnerRes = MakeRes("ASP", 5, new Vector3(5.0f, 0f, 0f));

        var graph = new InMemoryGraph();
        graph.Nodes[targetId]  = targetRes;
        graph.Nodes[partnerId] = partnerRes;

        // Electrostatic edge (LYS +1, ASP -1 → charge_product = -1)
        graph.Edges.Add(new Edge
        {
            FromId     = targetId,
            ToId       = partnerId,
            Type       = EdgeType.Electrostatic,
            DistanceA  = 5.0,
            EnergyKcal = -1.0
        });

        graph.BuildAdjacencyList();
        return (graph, targetId, partnerId);
    }

    /// <summary>
    /// Builds a 2-residue graph with a HYDROPHOBIC_CONTACT edge.
    /// Both residues are hydrophobic (ALA → target, LEU → partner).
    /// Used to test the BFS trace and the delta at hop 1.
    /// </summary>
    private static (InMemoryGraph graph, string targetId, string partnerId)
        BuildBfsGraph()
    {
        const string uid = "TEST";
        string targetId  = $"{uid}_A_1";
        string partnerId = $"{uid}_A_5";

        double mjAlaLeu = MiyazawaJernigan.GetEnergy("ALA", "LEU");

        var targetRes  = MakeRes("ALA", 1, new Vector3(0f,   0f, 0f), cb: new Vector3(0.5f,  0f, 0f));
        var partnerRes = MakeRes("LEU", 5, new Vector3(4.0f, 0f, 0f), cb: new Vector3(4.5f,  0f, 0f));

        var graph = new InMemoryGraph();
        graph.Nodes[targetId]  = targetRes;
        graph.Nodes[partnerId] = partnerRes;

        // Hydrophobic contact at CB-CB distance 4 Å
        graph.Edges.Add(new Edge
        {
            FromId     = targetId,
            ToId       = partnerId,
            Type       = EdgeType.Hydrophobic,
            DistanceA  = 4.0,
            EnergyKcal = mjAlaLeu
        });

        graph.BuildAdjacencyList();
        return (graph, targetId, partnerId);
    }

    // =========================================================================
    // 1. Residue node swap
    // =========================================================================

    [Fact]
    public void Swap_ResidueName_ChangedToMutantType()
    {
        var (graph, _, targetId, _) = BuildProGraph();
        var agent = new MutationAgent(graph);

        agent.ApplyMutationAndGetTrace(targetId, "PRO");

        Assert.Equal("PRO", graph.Nodes[targetId].ResidueName);
    }

    [Fact]
    public void Swap_OtherResidues_NotChanged()
    {
        var (graph, id1, targetId, id3) = BuildProGraph();
        var agent = new MutationAgent(graph);

        agent.ApplyMutationAndGetTrace(targetId, "PRO");

        Assert.Equal("ALA", graph.Nodes[id1].ResidueName);
        Assert.Equal("ASP", graph.Nodes[id3].ResidueName);
    }

    // =========================================================================
    // 2a. Edge rule: Any → Pro removes H-bond edges where target is donor
    // =========================================================================

    [Fact]
    public void ProMutation_RemovesHBondEdge_WhereMutatedNodeIsDonor()
    {
        var (graph, _, targetId, _) = BuildProGraph();
        var agent = new MutationAgent(graph);

        agent.ApplyMutationAndGetTrace(targetId, "PRO");

        bool donorHBondExists = graph.AdjacencyList[targetId]
            .Any(e => e.Type == EdgeType.HBond && e.Donor == targetId);

        Assert.False(donorHBondExists,
            "Mutating to PRO must remove H-bond edges where the mutated residue is the NH donor.");
    }

    [Fact]
    public void NonProMutation_PreservesHBondEdge()
    {
        var (graph, _, targetId, _) = BuildProGraph();
        var agent = new MutationAgent(graph);

        // Mutate to ALA (not PRO) — H-bond donor edge should be preserved
        agent.ApplyMutationAndGetTrace(targetId, "ALA");

        bool donorHBondExists = graph.AdjacencyList[targetId]
            .Any(e => e.Type == EdgeType.HBond && e.Donor == targetId);

        Assert.True(donorHBondExists,
            "Mutating to a non-PRO residue must leave H-bond donor edges intact.");
    }

    // =========================================================================
    // 2b. Edge rule: Any → hydrophobic creates new HYDROPHOBIC_CONTACT edges
    // =========================================================================

    [Fact]
    public void HydrophobicMutation_CreatesContact_WhenPartnerInRange()
    {
        // Target starts as GLU (not hydrophobic); partner is ILE (hydrophobic).
        // CB-CB distance = 4.0 Å < 8.0 Å — within contact range.
        var (graph, targetId, partnerId) = BuildHydrophobicGraph();
        var agent = new MutationAgent(graph);

        agent.ApplyMutationAndGetTrace(targetId, "LEU"); // LEU is hydrophobic

        bool contactExists = graph.AdjacencyList[targetId]
            .Any(e => e.Type == EdgeType.Hydrophobic
                   && (e.ToId == partnerId || e.FromId == partnerId));

        Assert.True(contactExists,
            "Mutating to a hydrophobic residue must create a HYDROPHOBIC_CONTACT edge " +
            "to a hydrophobic partner within 8 Å CB-CB range.");
    }

    [Fact]
    public void NonHydrophobicMutation_NoHydrophobicContact_Created()
    {
        // Target starts as GLU, mutate to ASP (also not hydrophobic).
        var (graph, targetId, partnerId) = BuildHydrophobicGraph();
        var agent = new MutationAgent(graph);

        agent.ApplyMutationAndGetTrace(targetId, "ASP");

        bool contactExists = graph.AdjacencyList[targetId]
            .Any(e => e.Type == EdgeType.Hydrophobic
                   && (e.ToId == partnerId || e.FromId == partnerId));

        Assert.False(contactExists,
            "Mutating to a non-hydrophobic residue must not create a HYDROPHOBIC_CONTACT edge.");
    }

    // =========================================================================
    // 2c. Edge rule: charged → neutral removes ELECTROSTATIC edges
    // =========================================================================

    [Fact]
    public void ChargeNeutralizingMutation_RemovesElectrostaticEdge()
    {
        // LYS (+1) has an electrostatic edge to ASP (-1); mutate LYS → ALA (0).
        var (graph, targetId, partnerId) = BuildElectrostaticGraph();
        var agent = new MutationAgent(graph);

        agent.ApplyMutationAndGetTrace(targetId, "ALA");

        bool electrostaticExists = graph.AdjacencyList[targetId]
            .Any(e => e.Type == EdgeType.Electrostatic
                   && (e.ToId == partnerId || e.FromId == partnerId));

        Assert.False(electrostaticExists,
            "Mutating a charged residue to a neutral one must remove the ELECTROSTATIC edge.");
    }

    // =========================================================================
    // 3. BFS energy propagation
    // =========================================================================

    [Fact]
    public void BfsTrace_IsNonEmpty_WhenGraphHasNeighbours()
    {
        var (graph, targetId, _) = BuildBfsGraph();
        var agent = new MutationAgent(graph);

        var (trace, _) = agent.ApplyMutationAndGetTrace(targetId, "VAL");

        Assert.NotEmpty(trace);
    }

    [Fact]
    public void BfsTrace_Hop1Partner_HasNonZeroDelta_WhenEnergyChanges()
    {
        // ALA → VAL changes the MJ energy with LEU partner.
        // The hop-1 trace entry for the partner must reflect this delta.
        var (graph, targetId, partnerId) = BuildBfsGraph();

        double expectedOld = MiyazawaJernigan.GetEnergy("ALA", "LEU");
        double expectedNew = MiyazawaJernigan.GetEnergy("VAL", "LEU");
        double expectedDelta = expectedNew - expectedOld;

        var agent = new MutationAgent(graph);
        var (trace, _) = agent.ApplyMutationAndGetTrace(targetId, "VAL");

        // The hop-1 trace entry for the partner captures the raw energy delta
        var hop1Entry = trace.FirstOrDefault(t => t.id == partnerId && t.hopDistance == 1);

        Assert.NotEqual(default, hop1Entry);
        Assert.NotEqual(0.0, hop1Entry.propagatedDelta);
        Assert.Equal(expectedDelta, hop1Entry.propagatedDelta, precision: 6);
    }

    [Fact]
    public void BfsTrace_IdenticalMutation_ProducesZeroDelta()
    {
        // Mutating ALA → ALA (no change) should give zero delta at the partner.
        var (graph, targetId, partnerId) = BuildBfsGraph();
        var agent = new MutationAgent(graph);

        var (trace, _) = agent.ApplyMutationAndGetTrace(targetId, "ALA");

        var hop1Entry = trace.FirstOrDefault(t => t.id == partnerId && t.hopDistance == 1);

        Assert.NotEqual(default, hop1Entry);
        Assert.Equal(0.0, hop1Entry.propagatedDelta, precision: 10);
    }
}

// =============================================================================
// 4. EnergySignalPropagator — interval update after each BFS step
// =============================================================================

public class EnergySignalPropagatorTests
{
    [Fact]
    public void Step1_LoAndHi_MatchExactFormula()
    {
        // After step 1: lo = DDG - ε₀, hi = DDG + ε₀  (because √1 = 1)
        double eps0  = 1.5;
        double delta = 0.8;
        var esp = new EnergySignalPropagator(eps0);

        esp.Apply(delta);

        Assert.Equal(delta - eps0, esp.Lo, precision: 12);
        Assert.Equal(delta + eps0, esp.Hi, precision: 12);
        Assert.Equal(2.0 * eps0,   esp.Width, precision: 12);
    }

    [Fact]
    public void MonotoneNarrowing_BetweenCrossovers_WidthNeverIncreases()
    {
        // Feed consistent positive deltas. lo can only increase (max), hi can only
        // decrease (min) → width must be non-increasing between crossover recovery events.
        // (A crossover recovery legitimately widens the interval to the unconstrained
        //  current estimate; the test skips width-increase checks for those steps.)
        var esp = new EnergySignalPropagator(epsilon0: 2.0);
        double lastWidth = double.PositiveInfinity;

        for (int k = 1; k <= 50; k++)
        {
            esp.Apply(0.1);

            if (!esp.CrossoverOccurred)
            {
                Assert.True(esp.Width <= lastWidth + 1e-12,
                    $"Width increased without crossover at step {k}: {lastWidth} → {esp.Width}");
            }

            lastWidth = esp.Width;
        }
    }

    [Fact]
    public void ZeroDeltas_Width_ConvergesMonotonicallyToZero()
    {
        // All deltas = 0.  DDG stays 0 throughout.
        // Uncertainty = eps0/√k → 0, so width = 2*eps0/√k → 0.
        var esp = new EnergySignalPropagator(epsilon0: 2.0);

        for (int k = 1; k <= 100; k++)
            esp.Apply(0.0);

        // After 100 steps: width = 2 * 2.0 / √100 = 0.4
        Assert.Equal(0.0, esp.CumulativeDDG, precision: 12);
        Assert.Equal(0.4, esp.Width, precision: 10);
        Assert.True(esp.Width < 4.0, "Width must have narrowed from initial 4.0 kcal/mol.");
    }

    [Fact]
    public void CumulativeDDG_IsCorrectRunningSum()
    {
        var esp = new EnergySignalPropagator(epsilon0: 1.0);
        double[] deltas = [0.5, -0.3, 1.2, -0.8, 0.4];
        double expected = 0;

        foreach (double d in deltas)
        {
            expected += d;
            esp.Apply(d);
        }

        Assert.Equal(expected, esp.CumulativeDDG, precision: 12);
        Assert.Equal(deltas.Length, esp.StepCount);
    }

    [Fact]
    public void CrossoverDetected_WhenLargeOscillatingDeltaCausesInversion()
    {
        // ε₀ = 0.05 (very tight). Step 1 locks [lo, hi] tightly around +2.
        // Step 2 swings to -5, which makes proposedHi < lo → crossover.
        var esp = new EnergySignalPropagator(epsilon0: 0.05);

        esp.Apply(2.0);   // step 1: lo ≈ 1.95, hi ≈ 2.05
        esp.Apply(-7.0);  // step 2: DDG = -5.0, proposedHi ≈ -4.96 < lo ≈ 1.95 → crossover

        Assert.True(esp.CrossoverOccurred,
            "A large opposing delta must trigger crossover recovery.");
    }

    [Fact]
    public void AfterCrossoverRecovery_IntervalIsValid()
    {
        // After crossover recovery the interval must still satisfy lo ≤ hi.
        var esp = new EnergySignalPropagator(epsilon0: 0.05);

        esp.Apply(2.0);
        esp.Apply(-7.0);

        Assert.True(esp.Lo <= esp.Hi,
            "After crossover recovery lo must be ≤ hi.");
    }

    [Fact]
    public void CrossoverFlag_ResetOnNextStep()
    {
        var esp = new EnergySignalPropagator(epsilon0: 0.05);
        esp.Apply(2.0);
        esp.Apply(-7.0);  // triggers crossover

        Assert.True(esp.CrossoverOccurred);

        esp.Apply(0.0);   // next step
        Assert.False(esp.CrossoverOccurred,
            "CrossoverOccurred must reset to false on the subsequent Apply() call.");
    }

    [Fact]
    public void InvalidEpsilon0_ThrowsArgumentOutOfRangeException()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new EnergySignalPropagator(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new EnergySignalPropagator(-1.5));
    }
}
