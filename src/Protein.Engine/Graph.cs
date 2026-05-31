using System.Collections.Generic;

namespace Protein.Engine;

public enum EdgeType
{
    Peptide,
    HBond,
    Disulfide,
    Hydrophobic,
    Electrostatic,
    VdW
}

public class Edge
{
    public string FromId { get; set; } = "";
    public string ToId { get; set; } = "";
    public EdgeType Type { get; set; }
    
    // Properties
    public double DistanceA { get; set; }
    public double EnergyKcal { get; set; }
    
    // Type specific
    public string? Donor { get; set; }
    public string? Acceptor { get; set; }
    public double DihedralDeg { get; set; }
}

public class InMemoryGraph
{
    public Dictionary<string, Residue> Nodes { get; } = new();
    public List<Edge> Edges { get; } = new();

    // Fast lookup
    public Dictionary<string, List<Edge>> AdjacencyList { get; } = new();

    public void BuildAdjacencyList()
    {
        AdjacencyList.Clear();
        foreach (var node in Nodes.Keys)
        {
            AdjacencyList[node] = new List<Edge>();
        }

        foreach (var edge in Edges)
        {
            if (AdjacencyList.TryGetValue(edge.FromId, out var listFrom))
                listFrom.Add(edge);
            
            if (AdjacencyList.TryGetValue(edge.ToId, out var listTo))
                listTo.Add(edge); // Undirected representation of incident edges
        }
    }

    /// <summary>
    /// Returns a deep clone of this graph so that mutation agents can modify
    /// their private copy without affecting the shared read-only original.
    /// </summary>
    public InMemoryGraph Clone()
    {
        var clone = new InMemoryGraph();

        // Clone residue nodes
        foreach (var (id, res) in Nodes)
        {
            clone.Nodes[id] = new Residue
            {
                SeqPos       = res.SeqPos,
                ResidueName  = res.ResidueName,
                ChainId      = res.ChainId,
                PLDDT        = res.PLDDT,
                N            = res.N,
                CA           = res.CA,
                C            = res.C,
                O            = res.O,
                CB           = res.CB,
                Phi          = res.Phi,
                Psi          = res.Psi,
                SideChainCentroid = res.SideChainCentroid
            };
        }

        // Clone edges
        foreach (var e in Edges)
        {
            clone.Edges.Add(new Edge
            {
                FromId     = e.FromId,
                ToId       = e.ToId,
                Type       = e.Type,
                DistanceA  = e.DistanceA,
                EnergyKcal = e.EnergyKcal,
                Donor      = e.Donor,
                Acceptor   = e.Acceptor,
                DihedralDeg = e.DihedralDeg
            });
        }

        clone.BuildAdjacencyList();
        return clone;
    }
}
