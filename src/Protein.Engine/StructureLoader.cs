using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Threading.Tasks;
using Neo4j.Driver;

namespace Protein.Engine;

public class StructureLoader : IAsyncDisposable
{
    private readonly IDriver _driver;

    public StructureLoader(string uri, string username, string password)
    {
        _driver = GraphDatabase.Driver(uri, AuthTokens.Basic(username, password));
    }

    /// <summary>
    /// Master orchestration method: reads .mmCIF, groups into residues, calculates dihedrals/centroids,
    /// detects edges, sets up the transaction-isolated database, and bulk-inserts nodes/edges.
    /// </summary>
    public async Task GraftFamilyAsync(string uniprotId, string familyName, string pfamId, string cifPath)
    {
        // 1. Parse mmCIF atoms
        Console.WriteLine($"Streaming atoms from {cifPath}...");
        var atoms = MmcifParser.StreamAtoms(cifPath).ToList();
        Console.WriteLine($"Parsed {atoms.Count} atom records.");

        if (atoms.Count == 0)
        {
            throw new InvalidOperationException($"No atoms parsed from mmCIF file: '{cifPath}'");
        }

        // 2. Group into residues
        var residues = BuildResidues(atoms);
        Console.WriteLine($"Grouped into {residues.Count} residues.");

        // 3. Compute dihedrals and centroids
        ComputeDihedralsAndCentroids(residues);

        // 4. Detect edges
        var edgePack = DetectEdges(residues, uniprotId);
        Console.WriteLine($"Detected edges:");
        Console.WriteLine($"  - Peptide Bonds (:PEPTIDE): {edgePack.Peptides.Count}");
        Console.WriteLine($"  - Hydrogen Bonds (:H_BOND): {edgePack.HBonds.Count}");
        Console.WriteLine($"  - Disulfide Bridges (:DISULFIDE): {edgePack.Disulfides.Count}");
        Console.WriteLine($"  - Hydrophobic Contacts (:HYDROPHOBIC_CONTACT): {edgePack.Hydrophobics.Count}");
        Console.WriteLine($"  - Electrostatic Pairs (:ELECTROSTATIC): {edgePack.Electrostatics.Count}");
        Console.WriteLine($"  - Van der Waals Interactions (:VAN_DER_WAALS): {edgePack.VdWs.Count}");

        // 5. Connect and setup database (dashed format, no underscores allowed)
        string databaseName = familyName.Replace("_", "-").ToLowerInvariant();
        Console.WriteLine($"Setting up Neo4j database: '{databaseName}'...");
        await SetupDatabaseAsync(databaseName);

        // 6. Bulk insert residues
        await InsertResiduesAsync(databaseName, residues, uniprotId, familyName, pfamId);

        // 7. Bulk insert edges
        await InsertEdgesAsync(databaseName, edgePack);

        // 8. Write audit record to run_registry
        await WriteGraftAuditAsync(
            familyName, uniprotId, pfamId,
            residues.Count,
            edgePack.Peptides.Count,
            edgePack.HBonds.Count,
            edgePack.Disulfides.Count,
            edgePack.Hydrophobics.Count,
            edgePack.Electrostatics.Count,
            edgePack.VdWs.Count);

        Console.WriteLine($"Successfully grafted family '{familyName}' into database '{databaseName}'!");
    }

    private List<Residue> BuildResidues(List<AtomSite> atoms)
    {
        var groups = atoms
            .GroupBy(a => new { a.ChainId, a.SeqPos })
            .OrderBy(g => g.Key.ChainId)
            .ThenBy(g => g.Key.SeqPos)
            .ToList();

        var residues = new List<Residue>();

        foreach (var g in groups)
        {
            var res = new Residue
            {
                SeqPos = g.Key.SeqPos,
                ResidueName = g.First().ResidueName,
                ChainId = g.Key.ChainId
            };

            res.AllAtoms.AddRange(g);

            foreach (var atom in g)
            {
                string name = atom.AtomName.Trim().ToUpperInvariant();
                var pos = new Vector3((float)atom.X, (float)atom.Y, (float)atom.Z);

                if (name == "N") res.N = pos;
                else if (name == "CA") res.CA = pos;
                else if (name == "C") res.C = pos;
                else if (name == "O") res.O = pos;
                else if (name == "CB") res.CB = pos;
            }

            // Set PLDDT to CA's TempFactor, fallback to first atom's TempFactor
            var caAtom = g.FirstOrDefault(a => a.AtomName.Trim().Equals("CA", StringComparison.OrdinalIgnoreCase));
            if (caAtom != null)
            {
                res.PLDDT = caAtom.TempFactor;
            }
            else
            {
                res.PLDDT = g.First().TempFactor;
            }

            residues.Add(res);
        }

        return residues;
    }

    private void ComputeDihedralsAndCentroids(List<Residue> residues)
    {
        var chains = residues.GroupBy(r => r.ChainId).ToList();

        foreach (var chain in chains)
        {
            var chainResidues = chain.OrderBy(r => r.SeqPos).ToList();
            int n = chainResidues.Count;

            for (int i = 0; i < n; i++)
            {
                var r = chainResidues[i];

                r.ComputeSideChainCentroid();

                // Calculate Phi (φ)
                if (i > 0)
                {
                    var rPrev = chainResidues[i - 1];
                    if (rPrev.C.HasValue && r.N.HasValue && r.CA.HasValue && r.C.HasValue)
                    {
                        double rad = CalculateDihedral(rPrev.C.Value, r.N.Value, r.CA.Value, r.C.Value);
                        r.Phi = rad * (180.0 / Math.PI);
                    }
                }

                // Calculate Psi (ψ)
                if (i < n - 1)
                {
                    var rNext = chainResidues[i + 1];
                    if (r.N.HasValue && r.CA.HasValue && r.C.HasValue && rNext.N.HasValue)
                    {
                        double rad = CalculateDihedral(r.N.Value, r.CA.Value, r.C.Value, rNext.N.Value);
                        r.Psi = rad * (180.0 / Math.PI);
                    }
                }
            }
        }
    }

    public static double CalculateDihedral(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 p4)
    {
        Vector3 b1 = p2 - p1;
        Vector3 b2 = p3 - p2;
        Vector3 b3 = p4 - p3;

        Vector3 n1 = Vector3.Cross(b1, b2);
        Vector3 n2 = Vector3.Cross(b2, b3);

        Vector3 b2Normalized = Vector3.Normalize(b2);
        Vector3 m1 = Vector3.Cross(n1, b2Normalized);

        double x = Vector3.Dot(m1, n2);
        double y = Vector3.Dot(n1, n2);

        return -Math.Atan2(x, y);
    }

    private EdgePack DetectEdges(List<Residue> residues, string uniprotId)
    {
        var pack = new EdgePack();
        int n = residues.Count;

        var interactionTracker = new HashSet<(string, string)>();

        void RegisterInteraction(string idA, string idB)
        {
            var pair = string.CompareOrdinal(idA, idB) < 0 ? (idA, idB) : (idB, idA);
            interactionTracker.Add(pair);
        }

        bool HasInteraction(string idA, string idB)
        {
            var pair = string.CompareOrdinal(idA, idB) < 0 ? (idA, idB) : (idB, idA);
            return interactionTracker.Contains(pair);
        }

        string GetId(Residue r) => $"{uniprotId}_{r.ChainId}_{r.SeqPos}";

        // 1. Peptide Bonds (consecutive in same chain)
        var chainGroups = residues.GroupBy(r => r.ChainId).ToList();
        foreach (var group in chainGroups)
        {
            var sorted = group.OrderBy(r => r.SeqPos).ToList();
            for (int i = 0; i < sorted.Count - 1; i++)
            {
                var r1 = sorted[i];
                var r2 = sorted[i + 1];

                if (r1.C.HasValue && r2.N.HasValue)
                {
                    double dist = Vector3.Distance(r1.C.Value, r2.N.Value);
                    if (dist <= 2.0)
                    {
                        string id1 = GetId(r1);
                        string id2 = GetId(r2);
                        pack.Peptides.Add(new PeptideEdge(id1, id2, dist));
                        RegisterInteraction(id1, id2);
                    }
                }
            }
        }

        // Pairwise checks
        for (int i = 0; i < n; i++)
        {
            var r1 = residues[i];
            string id1 = GetId(r1);

            for (int j = i + 1; j < n; j++)
            {
                var r2 = residues[j];
                string id2 = GetId(r2);

                bool isLinked = false;

                // 2. Disulfide Bridge (Both CYS, SG-SG distance <= 2.5)
                if (r1.ResidueName.Equals("CYS", StringComparison.OrdinalIgnoreCase) &&
                    r2.ResidueName.Equals("CYS", StringComparison.OrdinalIgnoreCase))
                {
                    var sg1 = r1.AllAtoms.FirstOrDefault(a => a.AtomName.Trim().Equals("SG", StringComparison.OrdinalIgnoreCase));
                    var sg2 = r2.AllAtoms.FirstOrDefault(a => a.AtomName.Trim().Equals("SG", StringComparison.OrdinalIgnoreCase));

                    if (sg1 != null && sg2 != null)
                    {
                        var pSg1 = new Vector3((float)sg1.X, (float)sg1.Y, (float)sg1.Z);
                        var pSg2 = new Vector3((float)sg2.X, (float)sg2.Y, (float)sg2.Z);
                        double distSg = Vector3.Distance(pSg1, pSg2);

                        if (distSg <= 2.5)
                        {
                            if (r1.CB.HasValue && r2.CB.HasValue)
                            {
                                double rad = CalculateDihedral(r1.CB.Value, pSg1, pSg2, r2.CB.Value);
                                double deg = rad * (180.0 / Math.PI);

                                double absDeg = Math.Abs(deg);
                                if (absDeg >= 60.0 && absDeg <= 120.0)
                                {
                                    double cbDist = Vector3.Distance(r1.CB.Value, r2.CB.Value);
                                    pack.Disulfides.Add(new DisulfideEdge(id1, id2, cbDist, deg));
                                    isLinked = true;
                                }
                            }
                        }
                    }
                }

                // 3. Hydrogen Bond (Any N/O/S on r1 to any N/O/S on r2 distance <= 3.5)
                var heavy1 = r1.AllAtoms.Where(a => a.TypeSymbol is "N" or "O" or "S" or "n" or "o" or "s").ToList();
                var heavy2 = r2.AllAtoms.Where(a => a.TypeSymbol is "N" or "O" or "S" or "n" or "o" or "s").ToList();

                double minHBondDist = double.MaxValue;
                AtomSite? bestA1 = null;
                AtomSite? bestA2 = null;

                foreach (var h1 in heavy1)
                {
                    var p1 = new Vector3((float)h1.X, (float)h1.Y, (float)h1.Z);
                    foreach (var h2 in heavy2)
                    {
                        var p2 = new Vector3((float)h2.X, (float)h2.Y, (float)h2.Z);
                        double d = Vector3.Distance(p1, p2);
                        if (d < minHBondDist)
                        {
                            minHBondDist = d;
                            bestA1 = h1;
                            bestA2 = h2;
                        }
                    }
                }

                if (minHBondDist <= 3.5 && bestA1 != null && bestA2 != null)
                {
                    string donor = id1;
                    string acceptor = id2;

                    if (bestA1.TypeSymbol.Equals("N", StringComparison.OrdinalIgnoreCase) && 
                        !bestA2.TypeSymbol.Equals("N", StringComparison.OrdinalIgnoreCase))
                    {
                        donor = id1;
                        acceptor = id2;
                    }
                    else if (bestA2.TypeSymbol.Equals("N", StringComparison.OrdinalIgnoreCase) && 
                             !bestA1.TypeSymbol.Equals("N", StringComparison.OrdinalIgnoreCase))
                    {
                        donor = id2;
                        acceptor = id1;
                    }

                    pack.HBonds.Add(new HBondEdge(id1, id2, -0.5, donor, acceptor));
                    isLinked = true;
                }

                // 4. Hydrophobic Contact (Both hydrophobic, CB-CB distance <= 8.0)
                if (MiyazawaJernigan.IsHydrophobic(r1.ResidueName) &&
                    MiyazawaJernigan.IsHydrophobic(r2.ResidueName))
                {
                    if (r1.CB.HasValue && r2.CB.HasValue)
                    {
                        double distCb = Vector3.Distance(r1.CB.Value, r2.CB.Value);
                        if (distCb <= 8.0)
                        {
                            double energy = MiyazawaJernigan.GetEnergy(r1.ResidueName, r2.ResidueName);
                            pack.Hydrophobics.Add(new HydrophobicEdge(id1, id2, energy));
                            isLinked = true;
                        }
                    }
                }

                // 5. Electrostatic Pair (At least one charged, sidechain centroid distance <= 10.0)
                double chg1 = MiyazawaJernigan.GetCharge(r1.ResidueName);
                double chg2 = MiyazawaJernigan.GetCharge(r2.ResidueName);

                if (chg1 != 0.0 || chg2 != 0.0)
                {
                    if (r1.SideChainCentroid.HasValue && r2.SideChainCentroid.HasValue)
                    {
                        double distCent = Vector3.Distance(r1.SideChainCentroid.Value, r2.SideChainCentroid.Value);
                        if (distCent <= 10.0)
                        {
                            pack.Electrostatics.Add(new ElectrostaticEdge(id1, id2, chg1 * chg2, distCent));
                            isLinked = true;
                        }
                    }
                }

                if (isLinked)
                {
                    RegisterInteraction(id1, id2);
                }
            }
        }

        // 6. Van der Waals (Any pair, CA-CA distance 3.0 to 8.0, and NOT already linked by other interaction types)
        for (int i = 0; i < n; i++)
        {
            var r1 = residues[i];
            string id1 = GetId(r1);

            for (int j = i + 1; j < n; j++)
            {
                var r2 = residues[j];
                string id2 = GetId(r2);

                if (!HasInteraction(id1, id2))
                {
                    if (r1.CA.HasValue && r2.CA.HasValue)
                    {
                        double distCa = Vector3.Distance(r1.CA.Value, r2.CA.Value);
                        if (distCa >= 3.0 && distCa <= 8.0)
                        {
                            double energy = MiyazawaJernigan.GetEnergy(r1.ResidueName, r2.ResidueName);
                            pack.VdWs.Add(new VdWEdge(id1, id2, energy));
                        }
                    }
                }
            }
        }

        return pack;
    }

    private async Task SetupDatabaseAsync(string databaseName)
    {
        await using var session = _driver.AsyncSession(o => o.WithDatabase("system"));
        try
        {
            await session.ExecuteWriteAsync(async tx =>
            {
                await tx.RunAsync($"CREATE DATABASE `{databaseName}` IF NOT EXISTS");
            });
            Console.WriteLine($"Sent CREATE DATABASE command for '{databaseName}'.");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Warning / Info: Error during CREATE DATABASE: {ex.Message}");
        }

        // Wait for database to be online
        bool isOnline = false;
        for (int attempt = 0; attempt < 15; attempt++)
        {
            try
            {
                await using var checkSession = _driver.AsyncSession(o => o.WithDatabase("system"));
                var result = await checkSession.ExecuteReadAsync(async tx =>
                {
                    var cursor = await tx.RunAsync("SHOW DATABASES");
                    while (await cursor.FetchAsync())
                    {
                        var name = cursor.Current["name"].As<string>();
                        if (name.Equals(databaseName, StringComparison.OrdinalIgnoreCase))
                        {
                            var status = cursor.Current["currentStatus"].As<string>();
                            return status.Equals("online", StringComparison.OrdinalIgnoreCase);
                        }
                    }
                    return false;
                });
                if (result)
                {
                    isOnline = true;
                    break;
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error checking database status (attempt {attempt + 1}): {ex.Message}");
            }
            await Task.Delay(1000);
        }

        if (!isOnline)
        {
            Console.WriteLine($"Warning: Database '{databaseName}' was not reported as online. We will attempt to connect anyway.");
        }
        else
        {
            Console.WriteLine($"Database '{databaseName}' is ONLINE.");
        }
    }

    private async Task InsertResiduesAsync(
        string databaseName,
        List<Residue> residues,
        string uniprotId,
        string familyName,
        string pfamId)
    {
        Console.WriteLine($"Bulk inserting {residues.Count} residues into database '{databaseName}'...");
        await using var session = _driver.AsyncSession(o => o.WithDatabase(databaseName));

        var residueParams = residues.Select(r => new Dictionary<string, object?>
        {
            { "id", $"{uniprotId}_{r.ChainId}_{r.SeqPos}" },
            { "seqPos", r.SeqPos },
            { "name", r.ResidueName },
            { "chainId", r.ChainId },
            { "phi", r.Phi },
            { "psi", r.Psi },
            { "plddt", r.PLDDT },
            { "uniprotId", uniprotId },
            { "familyName", familyName },
            { "pfamId", pfamId },
            // Backbone + Cβ coordinates — stored as flat doubles; null when absent (e.g. Cβ on GLY)
            { "n_x",  (object?)(r.N  is { } n  ? (double?)n.X  : null) },
            { "n_y",  (object?)(r.N  is { } n2 ? (double?)n2.Y : null) },
            { "n_z",  (object?)(r.N  is { } n3 ? (double?)n3.Z : null) },
            { "ca_x", (object?)(r.CA is { } ca  ? (double?)ca.X  : null) },
            { "ca_y", (object?)(r.CA is { } ca2 ? (double?)ca2.Y : null) },
            { "ca_z", (object?)(r.CA is { } ca3 ? (double?)ca3.Z : null) },
            { "cb_x", (object?)(r.CB is { } cb  ? (double?)cb.X  : null) },
            { "cb_y", (object?)(r.CB is { } cb2 ? (double?)cb2.Y : null) },
            { "cb_z", (object?)(r.CB is { } cb3 ? (double?)cb3.Z : null) },
            { "sc_x", (object?)(r.SideChainCentroid is { } sc  ? (double?)sc.X  : null) },
            { "sc_y", (object?)(r.SideChainCentroid is { } sc2 ? (double?)sc2.Y : null) },
            { "sc_z", (object?)(r.SideChainCentroid is { } sc3 ? (double?)sc3.Z : null) },
        }).ToList();

        string query = @"
            UNWIND $residues AS row
            MERGE (r:Residue { id: row.id })
            SET r.seqPos = row.seqPos,
                r.name = row.name,
                r.chainId = row.chainId,
                r.phi = row.phi,
                r.psi = row.psi,
                r.plddt = row.plddt,
                r.uniprotId = row.uniprotId,
                r.familyName = row.familyName,
                r.pfamId = row.pfamId,
                r.n_x = row.n_x, r.n_y = row.n_y, r.n_z = row.n_z,
                r.ca_x = row.ca_x, r.ca_y = row.ca_y, r.ca_z = row.ca_z,
                r.cb_x = row.cb_x, r.cb_y = row.cb_y, r.cb_z = row.cb_z,
                r.sc_x = row.sc_x, r.sc_y = row.sc_y, r.sc_z = row.sc_z
        ";

        await session.ExecuteWriteAsync(async tx =>
        {
            await tx.RunAsync(query, new { residues = residueParams });
        });

        // Add unique constraint for residue ID
        try
        {
            await session.ExecuteWriteAsync(async tx =>
            {
                await tx.RunAsync("CREATE CONSTRAINT residue_id_unique IF NOT EXISTS FOR (r:Residue) REQUIRE r.id IS UNIQUE");
            });
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Constraint creation info: {ex.Message}");
        }
    }

    private async Task InsertEdgesAsync(string databaseName, EdgePack pack)
    {
        await using var session = _driver.AsyncSession(o => o.WithDatabase(databaseName));

        // 1. Peptide Bonds
        if (pack.Peptides.Count > 0)
        {
            Console.WriteLine($"Inserting {pack.Peptides.Count} peptide bonds...");
            var list = pack.Peptides.Select(e => new Dictionary<string, object>
            {
                { "fromId", e.FromId },
                { "toId", e.ToId },
                { "length_A", e.LengthA }
            }).ToList();

            string query = @"
                UNWIND $edges AS row
                MATCH (a:Residue { id: row.fromId })
                MATCH (b:Residue { id: row.toId })
                MERGE (a)-[e:PEPTIDE]->(b)
                SET e.length_A = row.length_A
            ";
            await session.ExecuteWriteAsync(async tx => await tx.RunAsync(query, new { edges = list }));
        }

        // 2. Hydrogen Bonds
        if (pack.HBonds.Count > 0)
        {
            Console.WriteLine($"Inserting {pack.HBonds.Count} hydrogen bonds...");
            var list = pack.HBonds.Select(e => new Dictionary<string, object>
            {
                { "fromId", e.FromId },
                { "toId", e.ToId },
                { "energy_kcal", e.EnergyKcal },
                { "donor", e.Donor },
                { "acceptor", e.Acceptor }
            }).ToList();

            string query = @"
                UNWIND $edges AS row
                MATCH (a:Residue { id: row.fromId })
                MATCH (b:Residue { id: row.toId })
                MERGE (a)-[e:H_BOND]->(b)
                SET e.energy_kcal = row.energy_kcal,
                    e.donor = row.donor,
                    e.acceptor = row.acceptor
            ";
            await session.ExecuteWriteAsync(async tx => await tx.RunAsync(query, new { edges = list }));
        }

        // 3. Disulfide Bridges
        if (pack.Disulfides.Count > 0)
        {
            Console.WriteLine($"Inserting {pack.Disulfides.Count} disulfide bridges...");
            var list = pack.Disulfides.Select(e => new Dictionary<string, object>
            {
                { "fromId", e.FromId },
                { "toId", e.ToId },
                { "cbeta_dist_A", e.CbetaDistA },
                { "dihedral_deg", e.DihedralDeg }
            }).ToList();

            string query = @"
                UNWIND $edges AS row
                MATCH (a:Residue { id: row.fromId })
                MATCH (b:Residue { id: row.toId })
                MERGE (a)-[e:DISULFIDE]->(b)
                SET e.cbeta_dist_A = row.cbeta_dist_A,
                    e.dihedral_deg = row.dihedral_deg
            ";
            await session.ExecuteWriteAsync(async tx => await tx.RunAsync(query, new { edges = list }));
        }

        // 4. Hydrophobic Contacts
        if (pack.Hydrophobics.Count > 0)
        {
            Console.WriteLine($"Inserting {pack.Hydrophobics.Count} hydrophobic contacts...");
            var list = pack.Hydrophobics.Select(e => new Dictionary<string, object>
            {
                { "fromId", e.FromId },
                { "toId", e.ToId },
                { "mj_energy_kcal", e.MjEnergyKcal }
            }).ToList();

            string query = @"
                UNWIND $edges AS row
                MATCH (a:Residue { id: row.fromId })
                MATCH (b:Residue { id: row.toId })
                MERGE (a)-[e:HYDROPHOBIC_CONTACT]->(b)
                SET e.mj_energy_kcal = row.mj_energy_kcal
            ";
            await session.ExecuteWriteAsync(async tx => await tx.RunAsync(query, new { edges = list }));
        }

        // 5. Electrostatic Pairs
        if (pack.Electrostatics.Count > 0)
        {
            Console.WriteLine($"Inserting {pack.Electrostatics.Count} electrostatic pairs...");
            var list = pack.Electrostatics.Select(e => new Dictionary<string, object>
            {
                { "fromId", e.FromId },
                { "toId", e.ToId },
                { "charge_product", e.ChargeProduct },
                { "distance_A", e.DistanceA }
            }).ToList();

            string query = @"
                UNWIND $edges AS row
                MATCH (a:Residue { id: row.fromId })
                MATCH (b:Residue { id: row.toId })
                MERGE (a)-[e:ELECTROSTATIC]->(b)
                SET e.charge_product = row.charge_product,
                    e.distance_A = row.distance_A
            ";
            await session.ExecuteWriteAsync(async tx => await tx.RunAsync(query, new { edges = list }));
        }

        // 6. Van der Waals
        if (pack.VdWs.Count > 0)
        {
            Console.WriteLine($"Inserting {pack.VdWs.Count} Van der Waals interactions...");
            var list = pack.VdWs.Select(e => new Dictionary<string, object>
            {
                { "fromId", e.FromId },
                { "toId", e.ToId },
                { "mj_energy_kcal", e.MjEnergyKcal }
            }).ToList();

            string query = @"
                UNWIND $edges AS row
                MATCH (a:Residue { id: row.fromId })
                MATCH (b:Residue { id: row.toId })
                MERGE (a)-[e:VAN_DER_WAALS]->(b)
                SET e.mj_energy_kcal = row.mj_energy_kcal
            ";
            await session.ExecuteWriteAsync(async tx => await tx.RunAsync(query, new { edges = list }));
        }
    }

    private async Task WriteGraftAuditAsync(
        string familyName,
        string uniprotId,
        string pfamId,
        int residueCount,
        int peptideCount,
        int hbondCount,
        int disulfideCount,
        int hydrophobicCount,
        int electrostaticCount,
        int vdwCount)
    {
        Console.WriteLine($"Writing graft audit to run-registry for '{familyName}'...");
        await using var session = _driver.AsyncSession(o => o.WithDatabase("run-registry"));

        await session.ExecuteWriteAsync(async tx =>
        {
            await tx.RunAsync(@"
                MERGE (a:GraftAudit { familyName: $familyName })
                SET a.uniprotId          = $uniprotId,
                    a.pfamId             = $pfamId,
                    a.residueCount       = $residueCount,
                    a.peptideCount       = $peptideCount,
                    a.hbondCount         = $hbondCount,
                    a.disulfideCount     = $disulfideCount,
                    a.hydrophobicCount   = $hydrophobicCount,
                    a.electrostaticCount = $electrostaticCount,
                    a.vdwCount           = $vdwCount,
                    a.timestampUtc       = $timestamp,
                    a.status             = 'complete'
            ", new
            {
                familyName,
                uniprotId,
                pfamId,
                residueCount,
                peptideCount,
                hbondCount,
                disulfideCount,
                hydrophobicCount,
                electrostaticCount,
                vdwCount,
                timestamp = DateTimeOffset.UtcNow.ToString("o")
            });
        });

        Console.WriteLine($"Audit record written for '{familyName}'.");
    }

    public async ValueTask DisposeAsync()
    {
        await _driver.DisposeAsync();
    }

    private class EdgePack
    {
        public List<PeptideEdge> Peptides { get; } = new();
        public List<HBondEdge> HBonds { get; } = new();
        public List<DisulfideEdge> Disulfides { get; } = new();
        public List<HydrophobicEdge> Hydrophobics { get; } = new();
        public List<ElectrostaticEdge> Electrostatics { get; } = new();
        public List<VdWEdge> VdWs { get; } = new();
    }

    private record PeptideEdge(string FromId, string ToId, double LengthA);
    private record HBondEdge(string FromId, string ToId, double EnergyKcal, string Donor, string Acceptor);
    private record DisulfideEdge(string FromId, string ToId, double CbetaDistA, double DihedralDeg);
    private record HydrophobicEdge(string FromId, string ToId, double MjEnergyKcal);
    private record ElectrostaticEdge(string FromId, string ToId, double ChargeProduct, double DistanceA);
    private record VdWEdge(string FromId, string ToId, double MjEnergyKcal);
}
