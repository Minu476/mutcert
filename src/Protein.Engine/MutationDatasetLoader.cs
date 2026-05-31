using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using Neo4j.Driver;

namespace Protein.Engine;

public class MutationDatasetLoader : IAsyncDisposable
{
    private readonly IDriver _driver;

    public MutationDatasetLoader(string uri, string username, string password)
    {
        _driver = GraphDatabase.Driver(uri, AuthTokens.Basic(username, password));
    }

    private static readonly Dictionary<char, string> OneToThreeLetter = new()
    {
        {'A', "ALA"}, {'R', "ARG"}, {'N', "ASN"}, {'D', "ASP"},
        {'C', "CYS"}, {'Q', "GLN"}, {'E', "GLU"}, {'G', "GLY"},
        {'H', "HIS"}, {'I', "ILE"}, {'L', "LEU"}, {'K', "LYS"},
        {'M', "MET"}, {'F', "PHE"}, {'P', "PRO"}, {'S', "SER"},
        {'T', "THR"}, {'W', "TRP"}, {'Y', "TYR"}, {'V', "VAL"}
    };

    public async Task ImportMutationsAsync(string csvPath)
    {
        var targetFamilies = new Dictionary<string, string>
        {
            {"P00720", "t4-lysozyme"},
            {"P01053", "ci2"},
            {"P00648", "barnase"}
        };

        var mutationsByDb = new Dictionary<string, List<Dictionary<string, object>>>();
        foreach (var db in targetFamilies.Values)
        {
            mutationsByDb[db] = new List<Dictionary<string, object>>();
        }

        Console.WriteLine($"Reading mutations from {csvPath}...");
        
        using var reader = new StreamReader(csvPath);
        string? header = await reader.ReadLineAsync(); // skip header

        int lineCount = 0;
        int importedCount = 0;
        while (await reader.ReadLineAsync() is { } line)
        {
            lineCount++;
            var parts = line.Split(',');
            if (parts.Length < 7) continue;

            string id = parts[0].Trim();
            string uniprot = parts[1].Trim();
            string pdbWild = parts[2].Trim();
            string mutatedChain = parts[3].Trim();
            // skip parts[4] which is pos (normalized to 0)
            string mutationCode = parts[5].Trim();
            string ddgStr = parts[6].Trim();

            if (!targetFamilies.TryGetValue(uniprot, out string? dbName))
                continue;

            if (!double.TryParse(ddgStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double ddg))
                continue;

            // Parse mutationCode (e.g. "W138Y")
            var match = Regex.Match(mutationCode, @"^([A-Z])(\d+)([A-Z])$", RegexOptions.IgnoreCase);
            if (!match.Success)
            {
                Console.WriteLine($"Warning: Invalid mutation code '{mutationCode}' on line {lineCount + 1}");
                continue;
            }

            char wildChar = char.ToUpperInvariant(match.Groups[1].Value[0]);
            int seqPos = int.Parse(match.Groups[2].Value);
            char mutChar = char.ToUpperInvariant(match.Groups[3].Value[0]);

            if (!OneToThreeLetter.TryGetValue(wildChar, out string? wildType) ||
                !OneToThreeLetter.TryGetValue(mutChar, out string? mutantType))
            {
                Console.WriteLine($"Warning: Unknown amino acid in code '{mutationCode}' on line {lineCount + 1}");
                continue;
            }

            string mutationId = $"MUT_{uniprot}_{mutatedChain}_{seqPos}_{mutantType}";
            string residueId = $"{uniprot}_{mutatedChain}_{seqPos}";

            var param = new Dictionary<string, object>
            {
                { "id", mutationId },
                { "uniprot", uniprot },
                { "pdbWild", pdbWild },
                { "mutatedChain", mutatedChain },
                { "seqPos", seqPos },
                { "wildType", wildType },
                { "mutantType", mutantType },
                { "ddg", ddg },
                { "residueId", residueId }
            };

            mutationsByDb[dbName].Add(param);
            importedCount++;
        }

        Console.WriteLine($"Parsed {importedCount} target mutations across {targetFamilies.Count} families.");

        foreach (var (dbName, mutations) in mutationsByDb)
        {
            if (mutations.Count == 0) continue;

            Console.WriteLine($"Inserting {mutations.Count} mutations into '{dbName}'...");
            await using var session = _driver.AsyncSession(o => o.WithDatabase(dbName));
            
            string query = @"
                UNWIND $mutations AS row
                MERGE (m:Mutation { id: row.id })
                SET m.uniprot = row.uniprot,
                    m.pdbWild = row.pdbWild,
                    m.mutatedChain = row.mutatedChain,
                    m.seqPos = row.seqPos,
                    m.wildType = row.wildType,
                    m.mutantType = row.mutantType,
                    m.ddg_kcal_mol = row.ddg
                WITH m, row
                MATCH (r:Residue { id: row.residueId })
                MERGE (r)-[:HAS_MUTATION]->(m)
            ";

            await session.ExecuteWriteAsync(async tx =>
            {
                await tx.RunAsync(query, new { mutations });
            });
            
            Console.WriteLine($"Successfully loaded {mutations.Count} mutations into '{dbName}'.");
        }
    }

    public async ValueTask DisposeAsync()
    {
        await _driver.DisposeAsync();
    }
}
