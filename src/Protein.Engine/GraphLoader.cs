using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading.Tasks;
using Neo4j.Driver;

namespace Protein.Engine;

public static class GraphLoader
{
    private static readonly Dictionary<string, string> UniprotToDb = new()
    {
        {"P00720", "t4-lysozyme"},
        {"P01053", "ci2"},
        {"P00648", "barnase"}
    };

    public static async Task<InMemoryGraph> LoadGraphAsync(string uniprotId, string uri, string username, string password)
    {
        if (!UniprotToDb.TryGetValue(uniprotId, out string? dbName))
        {
            throw new ArgumentException($"Unknown UniProt ID: {uniprotId}");
        }

        var graph = new InMemoryGraph();
        await using var driver = GraphDatabase.Driver(uri, AuthTokens.Basic(username, password));
        await using var session = driver.AsyncSession(o => o.WithDatabase(dbName));

        // Load Nodes
        await session.ExecuteReadAsync(async tx =>
        {
            var cursor = await tx.RunAsync("MATCH (r:Residue {uniprotId: $u}) RETURN r", new { u = uniprotId });
            while (await cursor.FetchAsync())
            {
                var node = cursor.Current["r"].As<INode>();
                var props = node.Properties;

                string id = props["id"].As<string>();
                var res = new Residue
                {
                    SeqPos      = props["seqPos"].As<int>(),
                    ResidueName = props["name"].As<string>(),
                    ChainId     = props["chainId"].As<string>(),
                    PLDDT       = props["plddt"].As<double>(),
                    N           = ReadVec3(props, "n"),
                    CA          = ReadVec3(props, "ca"),
                    CB          = ReadVec3(props, "cb"),
                    SideChainCentroid = ReadVec3(props, "sc"),
                };

                graph.Nodes[id] = res;
            }
        });
        
        // Load Edges
        await session.ExecuteReadAsync(async tx =>
        {
            var cursor = await tx.RunAsync("MATCH (a:Residue {uniprotId: $u})-[e]->(b:Residue {uniprotId: $u}) RETURN a.id as fromId, type(e) as type, e, b.id as toId", new { u = uniprotId });
            while (await cursor.FetchAsync())
            {
                string fromId = cursor.Current["fromId"].As<string>();
                string toId = cursor.Current["toId"].As<string>();
                string typeStr = cursor.Current["type"].As<string>();
                var eProps = cursor.Current["e"].As<IRelationship>().Properties;
                
                EdgeType type = typeStr switch
                {
                    "PEPTIDE" => EdgeType.Peptide,
                    "H_BOND" => EdgeType.HBond,
                    "DISULFIDE" => EdgeType.Disulfide,
                    "HYDROPHOBIC_CONTACT" => EdgeType.Hydrophobic,
                    "ELECTROSTATIC" => EdgeType.Electrostatic,
                    "VAN_DER_WAALS" => EdgeType.VdW,
                    _ => throw new Exception($"Unknown edge type {typeStr}")
                };

                var edge = new Edge
                {
                    FromId = fromId,
                    ToId = toId,
                    Type = type,
                    DistanceA = eProps.ContainsKey("distance_A") ? eProps["distance_A"].As<double>() : 
                                eProps.ContainsKey("length_A") ? eProps["length_A"].As<double>() : 
                                eProps.ContainsKey("cbeta_dist_A") ? eProps["cbeta_dist_A"].As<double>() : 0.0,
                    EnergyKcal = eProps.ContainsKey("energy_kcal") ? eProps["energy_kcal"].As<double>() : 
                                 eProps.ContainsKey("mj_energy_kcal") ? eProps["mj_energy_kcal"].As<double>() : 0.0,
                    Donor = eProps.ContainsKey("donor") ? eProps["donor"].As<string>() : null,
                    Acceptor = eProps.ContainsKey("acceptor") ? eProps["acceptor"].As<string>() : null,
                    DihedralDeg = eProps.ContainsKey("dihedral_deg") ? eProps["dihedral_deg"].As<double>() : 0.0
                };
                
                graph.Edges.Add(edge);
            }
        });
        
        graph.BuildAdjacencyList();

        return graph;
    }

    /// <summary>
    /// Reads a Vector3 from Neo4j node properties stored as three separate double fields
    /// named <c>{prefix}_x</c>, <c>{prefix}_y</c>, <c>{prefix}_z</c>.
    /// Returns null if any component is missing or null in the database
    /// (e.g. Cβ is absent for glycine residues).
    /// </summary>
    private static Vector3? ReadVec3(
        IReadOnlyDictionary<string, object> props,
        string prefix)
    {
        string kx = $"{prefix}_x", ky = $"{prefix}_y", kz = $"{prefix}_z";
        if (!props.TryGetValue(kx, out var rx) || rx is null) return null;
        if (!props.TryGetValue(ky, out var ry) || ry is null) return null;
        if (!props.TryGetValue(kz, out var rz) || rz is null) return null;
        try
        {
            float x = (float)Convert.ToDouble(rx);
            float y = (float)Convert.ToDouble(ry);
            float z = (float)Convert.ToDouble(rz);
            return new Vector3(x, y, z);
        }
        catch
        {
            return null;
        }
    }
}
