using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Protein.Engine;

public class MutationAgent
{
    public InMemoryGraph Graph { get; }

    public MutationAgent(InMemoryGraph graph)
    {
        Graph = graph;
    }

    public (List<(string id, double propagatedDelta, int hopDistance)>, HashSet<string>) ApplyMutationAndGetTrace(string targetResidueId, string mutantType)
    {
        var rawDeltas = ComputeRawDeltasAndRewire(targetResidueId, mutantType);
        return RunPropagation(targetResidueId, rawDeltas);
    }

    private Dictionary<string, double> ComputeRawDeltasAndRewire(string targetResidueId, string mutantType)
    {
        var rawDeltas = new Dictionary<string, double>();
        if (!Graph.Nodes.TryGetValue(targetResidueId, out var targetNode))
            throw new ArgumentException($"Residue {targetResidueId} not found in graph.");

        string wildType = targetNode.ResidueName;
        targetNode.ResidueName = mutantType;

        // Recompute centroids if possible, though exact rotamers aren't built
        // We will just use the CA or CB for proxies since sidechain atoms haven't been re-sampled.
        
        // Edge Re-evaluations
        var targetEdges = Graph.AdjacencyList[targetResidueId].ToList();
        
        // Record original edges and their energies
        var originalEdges = targetEdges
            .GroupBy(e => e.FromId == targetResidueId ? e.ToId : e.FromId)
            .ToDictionary(g => g.Key, g => g.Sum(e => e.EnergyKcal));
        
        // 1. Any -> Pro: Delete hydrogen bonds where mutated residue is donor
        if (mutantType == "PRO")
        {
            var hBondsToRemove = targetEdges.Where(e => 
                e.Type == EdgeType.HBond && e.Donor == targetResidueId).ToList();
            
            foreach (var e in hBondsToRemove)
            {
                Graph.Edges.Remove(e);
                Graph.AdjacencyList[e.FromId].Remove(e);
                Graph.AdjacencyList[e.ToId].Remove(e);
            }
        }

        // 2. Any -> Cys / Cys -> Any: Re-evaluate Disulfide Bridges
        if (mutantType == "CYS" || wildType == "CYS")
        {
            if (mutantType != "CYS")
            {
                // Remove existing disulfides
                var ssToRemove = targetEdges.Where(e => e.Type == EdgeType.Disulfide).ToList();
                foreach (var e in ssToRemove)
                {
                    Graph.Edges.Remove(e);
                    Graph.AdjacencyList[e.FromId].Remove(e);
                    Graph.AdjacencyList[e.ToId].Remove(e);
                }
            }
            else
            {
                // Check for new CYS partners
                foreach (var (id, node) in Graph.Nodes)
                {
                    if (id != targetResidueId && node.ResidueName == "CYS")
                    {
                        if (targetNode.CB.HasValue && node.CB.HasValue)
                        {
                            double dist = Vector3.Distance(targetNode.CB.Value, node.CB.Value);
                            if (dist <= 7.5)
                            {
                                // Avoid duplicate
                                if (!Graph.AdjacencyList[targetResidueId].Any(e => 
                                    e.Type == EdgeType.Disulfide && (e.ToId == id || e.FromId == id)))
                                {
                                    var ssEdge = new Edge 
                                    {
                                        FromId = targetResidueId, ToId = id, Type = EdgeType.Disulfide, DistanceA = dist
                                    };
                                    Graph.Edges.Add(ssEdge);
                                    Graph.AdjacencyList[targetResidueId].Add(ssEdge);
                                    Graph.AdjacencyList[id].Add(ssEdge);
                                }
                            }
                        }
                    }
                }
            }
        }

        // 3. Any -> Charged: Re-evaluate electrostatics
        double mutCharge = MiyazawaJernigan.GetCharge(mutantType);
        double wildCharge = MiyazawaJernigan.GetCharge(wildType);
        if (mutCharge != wildCharge)
        {
            // Remove old electrostatics
            var elecToRemove = targetEdges.Where(e => e.Type == EdgeType.Electrostatic).ToList();
            foreach (var e in elecToRemove)
            {
                Graph.Edges.Remove(e);
                Graph.AdjacencyList[e.FromId].Remove(e);
                Graph.AdjacencyList[e.ToId].Remove(e);
            }

            if (mutCharge != 0.0)
            {
                // Add new electrostatics
                foreach (var (id, node) in Graph.Nodes)
                {
                    if (id == targetResidueId) continue;
                    double nCharge = MiyazawaJernigan.GetCharge(node.ResidueName);
                    if (nCharge != 0.0)
                    {
                        var p1 = targetNode.SideChainCentroid ?? targetNode.CA;
                        var p2 = node.SideChainCentroid ?? node.CA;
                        if (p1.HasValue && p2.HasValue)
                        {
                            double dist = Vector3.Distance(p1.Value, p2.Value);
                            if (dist <= 10.0)
                            {
                                var elecEdge = new Edge 
                                {
                                    FromId = targetResidueId, ToId = id, Type = EdgeType.Electrostatic,
                                    DistanceA = dist, EnergyKcal = mutCharge * nCharge // proxy for EnergyKcal charge_product
                                };
                                Graph.Edges.Add(elecEdge);
                                Graph.AdjacencyList[targetResidueId].Add(elecEdge);
                                Graph.AdjacencyList[id].Add(elecEdge);
                            }
                        }
                    }
                }
            }
        }

        // 4. Any -> Gly: Re-evaluate VdW
        if (mutantType == "GLY" || wildType == "GLY")
        {
            // Simple approach: re-calculate VdW energy for all existing VdW edges connected to target
            var vdwEdges = Graph.AdjacencyList[targetResidueId].Where(e => e.Type == EdgeType.VdW).ToList();
            foreach (var e in vdwEdges)
            {
                string partnerId = e.FromId == targetResidueId ? e.ToId : e.FromId;
                string partnerName = Graph.Nodes[partnerId].ResidueName;
                e.EnergyKcal = MiyazawaJernigan.GetEnergy(mutantType, partnerName);
            }
        }

        // 5. Any -> Any: Recreate hydrophobic contacts
        var hydroToRemove = targetEdges.Where(e => e.Type == EdgeType.Hydrophobic).ToList();
        foreach (var e in hydroToRemove)
        {
            Graph.Edges.Remove(e);
            Graph.AdjacencyList[e.FromId].Remove(e);
            Graph.AdjacencyList[e.ToId].Remove(e);
        }

        if (MiyazawaJernigan.IsHydrophobic(mutantType))
        {
            foreach (var (id, node) in Graph.Nodes)
            {
                if (id == targetResidueId) continue;
                if (MiyazawaJernigan.IsHydrophobic(node.ResidueName))
                {
                    if (targetNode.CB.HasValue && node.CB.HasValue)
                    {
                        double dist = Vector3.Distance(targetNode.CB.Value, node.CB.Value);
                        if (dist <= 8.0)
                        {
                            var he = new Edge 
                            {
                                FromId = targetResidueId, ToId = id, Type = EdgeType.Hydrophobic,
                                DistanceA = dist, EnergyKcal = MiyazawaJernigan.GetEnergy(mutantType, node.ResidueName)
                            };
                            Graph.Edges.Add(he);
                            Graph.AdjacencyList[targetResidueId].Add(he);
                            Graph.AdjacencyList[id].Add(he);
                        }
                    }
                }
            }
        }
        // Calculate differences for all neighbors
        var newEdges = Graph.AdjacencyList[targetResidueId];
        var newEdgeDict = newEdges
            .GroupBy(e => e.FromId == targetResidueId ? e.ToId : e.FromId)
            .ToDictionary(g => g.Key, g => g.Sum(e => e.EnergyKcal));

        var allNeighbors = new HashSet<string>(originalEdges.Keys.Concat(newEdgeDict.Keys));
        foreach (var n in allNeighbors)
        {
            double oldE = originalEdges.TryGetValue(n, out var oe) ? oe : 0.0;
            double newE = newEdgeDict.TryGetValue(n, out var ne) ? ne : 0.0;
            rawDeltas[n] = newE - oldE;
        }

        return rawDeltas;
    }

    private (List<(string id, double propagatedDelta, int hopDistance)>, HashSet<string>) RunPropagation(string targetResidueId, Dictionary<string, double> rawDeltas)
    {
        var trace = new List<(string id, double propagatedDelta, int hopDistance)>();
        var shellHops = new HashSet<string>();
        
        if (!Graph.Nodes.TryGetValue(targetResidueId, out var targetNode)) 
            return (trace, shellHops);

        var visited = new HashSet<string>();
        var queue = new Queue<(string id, double propagatedDelta, int hopDistance)>();
        
        visited.Add(targetResidueId);
        
        // Find 2-hop neighborhood to represent the unconditionally evaluated structural shell
        var hop1 = Graph.AdjacencyList[targetResidueId].Select(e => e.FromId == targetResidueId ? e.ToId : e.FromId).ToHashSet();
        var hop2 = new HashSet<string>();
        foreach (var h1 in hop1)
        {
            foreach (var e in Graph.AdjacencyList[h1])
            {
                string next = e.FromId == h1 ? e.ToId : e.FromId;
                if (next != targetResidueId && !hop1.Contains(next)) hop2.Add(next);
            }
        }
        
        shellHops.UnionWith(hop1);
        shellHops.UnionWith(hop2);

        foreach (var (neighborId, delta) in rawDeltas)
        {
            queue.Enqueue((neighborId, delta, 1));
            visited.Add(neighborId);
            shellHops.Remove(neighborId);
        }

        // Enqueue remaining shell residues with delta = 0
        foreach (var id in shellHops)
        {
            queue.Enqueue((id, 0.0, hop1.Contains(id) ? 1 : 2));
            visited.Add(id);
        }
        
        // Add target node to shellHops originally just to have complete set tracking
        shellHops.UnionWith(hop1);
        shellHops.UnionWith(hop2);
        
        while (queue.Count > 0)
        {
            var (currentId, currentDelta, hopDist) = queue.Dequeue();
            trace.Add((currentId, currentDelta, hopDist));

            if (!Graph.Nodes.TryGetValue(currentId, out var currentNode)) continue;

            // Propagate to neighbors
            foreach (var edge in Graph.AdjacencyList[currentId])
            {
                string nextId = edge.FromId == currentId ? edge.ToId : edge.FromId;
                if (!visited.Contains(nextId))
                {
                    int nextHop = hopDist + 1;
                    double effectiveDelta = currentDelta * Math.Exp(-nextHop / 2.5);

                    if (nextHop <= 2 || Math.Abs(effectiveDelta) > 0.01)
                    {
                        visited.Add(nextId);
                        queue.Enqueue((nextId, effectiveDelta, nextHop));
                    }
                }
            }
        }
        
        return (trace, shellHops);
    }
}
