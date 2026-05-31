using System;
using System.Diagnostics;
using System.Threading.Tasks;
using Neo4j.Driver;

namespace Protein.Engine;

internal class Program
{
    private static async Task<int> Main(string[] args)
    {
        if (args.Length == 0)
        {
            PrintUsage();
            return 1;
        }

        string subcommand = args[0].ToLowerInvariant();

        if (subcommand == "verify-isolation")
        {
            return await RunVerifyIsolationAsync();
        }

        if (subcommand == "verify")
        {
            return await RunVerifyAsync();
        }

        if (subcommand == "import-mutations")
        {
            return await RunImportMutationsAsync();
        }

        if (subcommand == "run-mutation")
        {
            return await RunMutationAsync(args);
        }

        if (subcommand == "batch-run")
        {
            return await RunBatchAsync(args);
        }

        if (subcommand == "whereami")
        {
            return await RunWhereAmIAsync();
        }

        if (subcommand == "replay")
        {
            return await RunReplayAsync(args);
        }

        if (subcommand == "validate-family")
        {
            return await RunValidateFamilyAsync(args);
        }

        if (subcommand == "validate-all")
        {
            return await RunValidateAllAsync();
        }

        if (subcommand == "zero-forgetting")
        {
            return await RunZeroForgettingAsync();
        }

        if (subcommand == "trace-mutation")
        {
            return await RunTraceMutationAsync(args);
        }

        if (subcommand == "graft")
        {
            if (args.Length < 5)
            {
                Console.WriteLine("Error: Missing required arguments for 'graft' command.");
                PrintUsage();
                return 1;
            }

            string uniprotId = args[1];
            string familyName = args[2];
            string pfamId = args[3];
            string cifPath = args[4];

            string uri = Environment.GetEnvironmentVariable("NEO4J_URI") ?? "bolt://localhost:7687";
            string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
            string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

            Console.WriteLine("================================================================================");
            Console.WriteLine(" MutCert Phase 3 (Graph Construction) — Graft Command");
            Console.WriteLine("================================================================================");
            Console.WriteLine($"UniProt ID:  {uniprotId}");
            Console.WriteLine($"Family Name: {familyName}");
            Console.WriteLine($"Pfam ID:     {pfamId}");
            Console.WriteLine($"CIF Path:    {cifPath}");
            Console.WriteLine($"Neo4j URI:   {uri}");
            Console.WriteLine($"Username:    {username}");
            Console.WriteLine("================================================================================");

            var stopwatch = Stopwatch.StartNew();
            try
            {
                await using var loader = new StructureLoader(uri, username, password);
                await loader.GraftFamilyAsync(uniprotId, familyName, pfamId, cifPath);
                
                stopwatch.Stop();
                Console.WriteLine($"Graft completed successfully in {stopwatch.Elapsed.TotalSeconds:F2} seconds.");
                return 0;
            }
            catch (Exception ex)
            {
                stopwatch.Stop();
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine($"\nError during graft execution: {ex.Message}");
                Console.WriteLine(ex.StackTrace);
                Console.ResetColor();
                return 2;
            }
        }

        PrintUsage();
        return 1;
    }

    private static async Task<int> RunVerifyAsync()
    {
        string uri = Environment.GetEnvironmentVariable("NEO4J_URI") ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

        Console.WriteLine("================================================================================");
        Console.WriteLine(" MutCert Phase 3 (Graph Construction) — Verification Audit");
        Console.WriteLine("================================================================================");
        Console.WriteLine($"Neo4j URI: {uri}");
        Console.WriteLine($"Username:  {username}");
        Console.WriteLine("================================================================================");

        string[] dbs = ["t4-lysozyme", "ci2", "barnase"];
        
        try
        {
            using var driver = GraphDatabase.Driver(uri, AuthTokens.Basic(username, password));

            foreach (var db in dbs)
            {
                Console.WriteLine($"\nDatabase: {db}");
                Console.WriteLine("--------------------------------------------------------------------------------");
                try
                {
                    await using var session = driver.AsyncSession(o => o.WithDatabase(db));
                    
                    long residues = await GetCountAsync(session, "MATCH (r:Residue) RETURN count(r)");
                    long peptides = await GetCountAsync(session, "MATCH ()-[r:PEPTIDE]->() RETURN count(r)");
                    long hbonds = await GetCountAsync(session, "MATCH ()-[r:H_BOND]->() RETURN count(r)");
                    long disulfides = await GetCountAsync(session, "MATCH ()-[r:DISULFIDE]->() RETURN count(r)");
                    long hydrophobics = await GetCountAsync(session, "MATCH ()-[r:HYDROPHOBIC_CONTACT]->() RETURN count(r)");
                    long electrostatics = await GetCountAsync(session, "MATCH ()-[r:ELECTROSTATIC]->() RETURN count(r)");
                    long vdws = await GetCountAsync(session, "MATCH ()-[r:VAN_DER_WAALS]->() RETURN count(r)");

                    Console.WriteLine($"  Residues (nodes):        {residues}");
                    Console.WriteLine($"  Peptide Bonds:           {peptides}");
                    Console.WriteLine($"  Hydrogen Bonds:          {hbonds}");
                    Console.WriteLine($"  Disulfide Bridges:       {disulfides}");
                    Console.WriteLine($"  Hydrophobic Contacts:    {hydrophobics}");
                    Console.WriteLine($"  Electrostatic Pairs:     {electrostatics}");
                    Console.WriteLine($"  Van der Waals:           {vdws}");
                    Console.WriteLine($"  Total Connections:       {peptides + hbonds + disulfides + hydrophobics + electrostatics + vdws}");
                }
                catch (Exception ex)
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine($"  Error reading database: {ex.Message}");
                    Console.ResetColor();
                }
            }
        }
        catch (Exception ex)
        {
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"Fatal verification error: {ex.Message}");
            Console.ResetColor();
            return 2;
        }

        Console.WriteLine("\n================================================================================");
        return 0;
    }

    private static async Task<int> RunImportMutationsAsync()
    {
        string uri = Environment.GetEnvironmentVariable("NEO4J_URI") ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";
        string csvPath = "data/s2648/s2648.csv";

        Console.WriteLine("================================================================================");
        Console.WriteLine(" MutCert Milestone 1 — Import Mutations");
        Console.WriteLine("================================================================================");
        
        var stopwatch = Stopwatch.StartNew();
        try
        {
            await using var loader = new MutationDatasetLoader(uri, username, password);
            await loader.ImportMutationsAsync(csvPath);
            
            stopwatch.Stop();
            Console.WriteLine($"Import completed successfully in {stopwatch.Elapsed.TotalSeconds:F2} seconds.");
            return 0;
        }
        catch (Exception ex)
        {
            stopwatch.Stop();
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"\nError during mutation import: {ex.Message}");
            Console.WriteLine(ex.StackTrace);
            Console.ResetColor();
            return 2;
        }
    }

    private static async Task<int> RunMutationAsync(string[] args)
    {
        if (args.Length < 4)
        {
            Console.WriteLine("Error: Missing required arguments for 'run-mutation'.");
            Console.WriteLine("Usage: dotnet run -- run-mutation <uniprotId> <seqPos> <mutant3Letter>");
            return 1;
        }

        string uniprotId = args[1];
        int seqPos = int.Parse(args[2]);
        string mutantType = args[3].ToUpperInvariant();
        string chainId = "A"; // Assume chain A for simplistic tests
        string targetResidueId = $"{uniprotId}_{chainId}_{seqPos}";
        string mutationId = $"MUT_{uniprotId}_{chainId}_{seqPos}_{mutantType}";

        string uri = Environment.GetEnvironmentVariable("NEO4J_URI") ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

        Console.WriteLine("================================================================================");
        Console.WriteLine($" MutCert Milestone 4 — Run Mutation {mutationId}");
        Console.WriteLine("================================================================================");

        var stopwatch = Stopwatch.StartNew();
        try
        {
            Console.WriteLine($"Loading graph for {uniprotId}...");
            var graph = await GraphLoader.LoadGraphAsync(uniprotId, uri, username, password);
            Console.WriteLine($"Graph loaded: {graph.Nodes.Count} residues, {graph.Edges.Count} edges.");

            var agent = new MutationAgent(graph);
            
            Console.WriteLine("Executing BFS Propagation...");
            var (trace, shellHops) = agent.ApplyMutationAndGetTrace(targetResidueId, mutantType);

            // Calibrate ε₀ from the S2648 training data for this family.
            double epsilon0 = EpsilonCalibrator.Calibrate("data/s2648/s2648.csv", uniprotId);
            
            Console.WriteLine($"Passing {trace.Count} steps to ConvergenceSupervisor (ε₀ = {epsilon0:F4})...");
            var supervisor = new ConvergenceSupervisor(uri, username, password, epsilon0);
            var summary = await supervisor.ProcessRunAsync(mutationId, targetResidueId, shellHops, trace);

            Console.WriteLine($"Result: Converged={summary.Converged} | DDG={summary.FinalDDG:F3} [{summary.FinalLo:F3}, {summary.FinalHi:F3}] kcal/mol | Steps={summary.TotalSteps}");

            stopwatch.Stop();
            Console.WriteLine($"Run completed successfully in {stopwatch.Elapsed.TotalSeconds:F2} seconds.");
            return 0;
        }
        catch (Exception ex)
        {
            stopwatch.Stop();
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"\nError during run-mutation: {ex.Message}");
            Console.WriteLine(ex.StackTrace);
            Console.ResetColor();
            return 2;
        }
    }

    private static async Task<int> RunValidateFamilyAsync(string[] args)
    {
        if (args.Length < 2)
        {
            Console.WriteLine("Error: Missing required argument <familyName>.");
            return 1;
        }

        string familyName = args[1];
        string uri = Environment.GetEnvironmentVariable("NEO4J_URI") ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

        Console.WriteLine("================================================================================");
        Console.WriteLine($" MutCert Milestone 5 — Validate Family '{familyName}'");
        Console.WriteLine("================================================================================");

        var stopwatch = Stopwatch.StartNew();
        try
        {
            var pipeline = new ValidationPipeline(uri, username, password);
            await pipeline.RunValidationAsync(familyName);

            stopwatch.Stop();
            Console.WriteLine($"\nValidation completed successfully in {stopwatch.Elapsed.TotalSeconds:F2} seconds.");
            return 0;
        }
        catch (Exception ex)
        {
            stopwatch.Stop();
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"\nError during validation: {ex.Message}");
            Console.WriteLine(ex.StackTrace);
            Console.ResetColor();
            return 2;
        }
    }

    private static async Task<long> GetCountAsync(IAsyncSession session, string query)
    {
        return await session.ExecuteReadAsync(async tx =>
        {
            var cursor = await tx.RunAsync(query);
            if (await cursor.FetchAsync())
            {
                return cursor.Current[0].As<long>();
            }
            return 0L;
        });
    }

    private static async Task<int> RunVerifyIsolationAsync()
    {
        string uri = Environment.GetEnvironmentVariable("NEO4J_URI") ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

        Console.WriteLine("================================================================================");
        Console.WriteLine(" MutCert Phase 3 — Cross-Family Isolation Audit (run_registry)");
        Console.WriteLine("================================================================================");

        var familyMap = new Dictionary<string, string>
        {
            ["t4_lysozyme"] = "t4-lysozyme",
            ["ci2"]         = "ci2",
            ["barnase"]     = "barnase"
        };

        bool allPassed = true;

        try
        {
            using var driver = GraphDatabase.Driver(uri, AuthTokens.Basic(username, password));

            // 1. Read all GraftAudit entries from run_registry
            var audits = new Dictionary<string, (long res, long pep, long hb, long ds, long hy, long el, long vdw)>();

            await using var regSession = driver.AsyncSession(o => o.WithDatabase("run-registry"));
            await regSession.ExecuteReadAsync(async tx =>
            {
                var cursor = await tx.RunAsync(@"
                    MATCH (a:GraftAudit)
                    RETURN a.familyName   AS family,
                           a.residueCount AS res,
                           a.peptideCount AS pep,
                           a.hbondCount   AS hb,
                           a.disulfideCount AS ds,
                           a.hydrophobicCount AS hy,
                           a.electrostaticCount AS el,
                           a.vdwCount     AS vdw
                ");
                while (await cursor.FetchAsync())
                {
                    string fn = cursor.Current["family"].As<string>();
                    audits[fn] = (
                        cursor.Current["res"].As<long>(),
                        cursor.Current["pep"].As<long>(),
                        cursor.Current["hb"].As<long>(),
                        cursor.Current["ds"].As<long>(),
                        cursor.Current["hy"].As<long>(),
                        cursor.Current["el"].As<long>(),
                        cursor.Current["vdw"].As<long>()
                    );
                }
                return true;
            });

            if (audits.Count == 0)
            {
                Console.ForegroundColor = ConsoleColor.Yellow;
                Console.WriteLine("No GraftAudit entries found in run_registry.");
                Console.WriteLine("Run 'graft' commands first to generate audit records.");
                Console.ResetColor();
                return 1;
            }

            Console.WriteLine($"Found {audits.Count} audit entries in run_registry.\n");

            // 2. For each family, compare audit snapshot to live counts
            foreach (var (familyName, dbName) in familyMap)
            {
                if (!audits.TryGetValue(familyName, out var expected))
                {
                    Console.ForegroundColor = ConsoleColor.Yellow;
                    Console.WriteLine($"[SKIP] No audit entry found for family '{familyName}'.");
                    Console.ResetColor();
                    continue;
                }

                Console.WriteLine($"Family: {familyName}  →  database: {dbName}");

                await using var famSession = driver.AsyncSession(o => o.WithDatabase(dbName));
                long liveRes = await GetCountAsync(famSession, "MATCH (r:Residue) RETURN count(r)");
                long livePep = await GetCountAsync(famSession, "MATCH ()-[r:PEPTIDE]->() RETURN count(r)");
                long liveHb  = await GetCountAsync(famSession, "MATCH ()-[r:H_BOND]->() RETURN count(r)");
                long liveDs  = await GetCountAsync(famSession, "MATCH ()-[r:DISULFIDE]->() RETURN count(r)");
                long liveHy  = await GetCountAsync(famSession, "MATCH ()-[r:HYDROPHOBIC_CONTACT]->() RETURN count(r)");
                long liveEl  = await GetCountAsync(famSession, "MATCH ()-[r:ELECTROSTATIC]->() RETURN count(r)");
                long liveVdw = await GetCountAsync(famSession, "MATCH ()-[r:VAN_DER_WAALS]->() RETURN count(r)");

                bool match =
                    liveRes == expected.res &&
                    livePep == expected.pep &&
                    liveHb  == expected.hb  &&
                    liveDs  == expected.ds  &&
                    liveHy  == expected.hy  &&
                    liveEl  == expected.el  &&
                    liveVdw == expected.vdw;

                long liveExtraNodes = liveRes - expected.res;
                long liveExtraEdges = (livePep + liveHb + liveDs + liveHy + liveEl + liveVdw)
                                    - (expected.pep + expected.hb + expected.ds + expected.hy + expected.el + expected.vdw);

                Console.WriteLine($"  {"Item",-24} {"Audit":>8} {"Live":>8} {"Delta":>8}");
                Console.WriteLine($"  {"─────────────────────────────────────────────────────────────────────────────────"}");
                void Row(string label, long exp, long live) =>
                    Console.WriteLine($"  {label,-24} {exp,8} {live,8} {(live - exp >= 0 ? "+" : "")}{live - exp,7}");

                Row("Residues",            expected.res, liveRes);
                Row("Peptide Bonds",       expected.pep, livePep);
                Row("Hydrogen Bonds",      expected.hb,  liveHb);
                Row("Disulfide Bridges",   expected.ds,  liveDs);
                Row("Hydrophobic Contacts",expected.hy,  liveHy);
                Row("Electrostatic Pairs", expected.el,  liveEl);
                Row("Van der Waals",       expected.vdw, liveVdw);

                if (match)
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine($"  → PASS  Zero cross-family writes detected.");
                    Console.ResetColor();
                }
                else
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine($"  → FAIL  Unexpected delta: +{liveExtraNodes} nodes, +{liveExtraEdges} edges.");
                    Console.ResetColor();
                    allPassed = false;
                }
                Console.WriteLine();
            }
        }
        catch (Exception ex)
        {
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"Fatal isolation audit error: {ex.Message}");
            Console.ResetColor();
            return 2;
        }

        Console.WriteLine("================================================================================");
        if (allPassed)
        {
            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine(" ALL FAMILIES PASS — cross-family write isolation confirmed.");
            Console.ResetColor();
            return 0;
        }
        else
        {
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine(" ISOLATION VIOLATIONS DETECTED — see failures above.");
            Console.ResetColor();
            return 1;
        }
    }

    /// <summary>
    /// validate-all: Phase 9 — run validation across all 3 families, print
    /// calibration curves, and emit data/validation_report.html.
    /// </summary>
    private static async Task<int> RunValidateAllAsync()
    {
        string uri      = Environment.GetEnvironmentVariable("NEO4J_URI")      ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

        Console.WriteLine("================================================================================");
        Console.WriteLine(" MutCert Phase 9 — Validation (all families)");
        Console.WriteLine("================================================================================");

        var sw = Stopwatch.StartNew();
        try
        {
            var pipeline = new ValidationPipeline(uri, username, password);
            var results  = await pipeline.RunAllFamiliesAsync(
                csvPath:   "data/s2648/s2648.csv",
                splitFile: "data/s2648_split.json");

            foreach (var r in results)
                ValidationPipeline.PrintFamilyResult(r);

            string reportPath = "data/validation_report.html";
            ValidationReport.Generate(results, reportPath);

            sw.Stop();
            Console.WriteLine($"\nValidation completed in {sw.Elapsed.TotalMinutes:F1} min.");

            bool targetMet = results.Any(r =>
                r.CalibrationCurve.FirstOrDefault(c => c.WidthThreshold == 2.0) is { } row &&
                row.MutCertCoverage >= 0.80);

            if (targetMet)
            {
                Console.ForegroundColor = ConsoleColor.Green;
                Console.WriteLine("\n✓ PHASE 9 TARGET MET: >= 80% coverage at +-2.0 kcal/mol on at least 1 family.");
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Yellow;
                Console.WriteLine("\no Phase 9 target not yet met (< 80% coverage). See report for details.");
            }
            Console.ResetColor();
            return 0;
        }
        catch (Exception ex)
        {
            sw.Stop();
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"\nError during validate-all: {ex.Message}");
            Console.WriteLine(ex.StackTrace);
            Console.ResetColor();
            return 2;
        }
    }

    /// <summary>
    /// zero-forgetting: Phase 10 — sequentially graft all 3 families and verify
    /// that each graft leaves all other family databases unchanged.
    /// Saves data/zero_forgetting_report.html.
    /// </summary>
    private static async Task<int> RunZeroForgettingAsync()
    {
        string uri      = Environment.GetEnvironmentVariable("NEO4J_URI")      ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

        try
        {
            var result = await ZeroForgettingVerifier.RunAndReportAsync(
                uri, username, password,
                reportPath: "data/zero_forgetting_report.html");

            return result.AllPassed ? 0 : 1;
        }
        catch (Exception ex)
        {
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"\nFatal error during zero-forgetting: {ex.Message}");
            Console.WriteLine(ex.StackTrace);
            Console.ResetColor();
            return 2;
        }
    }

    // =========================================================================
    // Phase 11 — Glass-box trace
    // =========================================================================

    /// <summary>
    /// trace-mutation: loads a stored MutationRun from run-registry and
    /// generates a self-contained HTML glass-box causal-chain trace.
    /// Usage: dotnet run -- trace-mutation &lt;mutationId&gt; [outputPath]
    /// </summary>
    private static async Task<int> RunTraceMutationAsync(string[] args)
    {
        if (args.Length < 2)
        {
            Console.WriteLine("Error: trace-mutation requires <mutationId> [outputPath].");
            Console.WriteLine("Example: dotnet run -- trace-mutation MUT_P00720_A_102_LEU");
            return 1;
        }

        string mutationId  = args[1];
        string safeName    = mutationId.Replace('/', '_').Replace('\\', '_');
        string defaultOut  = Path.Combine("data", $"trace_{safeName}.html");
        string outputPath  = args.Length >= 3 ? args[2] : defaultOut;

        string uri      = Environment.GetEnvironmentVariable("NEO4J_URI")      ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

        Console.WriteLine("================================================================================");
        Console.WriteLine(" MutCert Phase 11 — Glass-Box Causal Trace");
        Console.WriteLine("================================================================================");
        Console.WriteLine($" Mutation ID : {mutationId}");
        Console.WriteLine($" Output file : {outputPath}");
        Console.WriteLine($" Neo4j URI   : {uri}");
        Console.WriteLine("================================================================================");

        var sw = Stopwatch.StartNew();
        try
        {
            await using var tracer = new GlassBoxTrace(uri, username, password);
            string outFile = await tracer.GenerateAsync(mutationId, outputPath);
            sw.Stop();
            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine($"\nTrace generated in {sw.Elapsed.TotalSeconds:F2}s → {outFile}");
            Console.ResetColor();
            return 0;
        }
        catch (Exception ex)
        {
            sw.Stop();
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"\nError: {ex.Message}");
            Console.ResetColor();
            return 2;
        }
    }

    private static void PrintUsage()
    {
        Console.WriteLine("MutCert CLI Tool — Graph Engine");
        Console.WriteLine();
        Console.WriteLine("Usage:");
        Console.WriteLine("  dotnet run -- graft <uniprotId> <familyName> <pfamId> <cifPath>");
        Console.WriteLine("  dotnet run -- verify");
        Console.WriteLine("  dotnet run -- verify-isolation");
        Console.WriteLine("  dotnet run -- import-mutations");
        Console.WriteLine("  dotnet run -- run-mutation <uniprotId> <seqPos> <mutant3Letter>");
        Console.WriteLine("  dotnet run -- batch-run <uniprotId> <familyDbName> [workerCount]");
        Console.WriteLine("  dotnet run -- whereami");
        Console.WriteLine("  dotnet run -- replay <runId>");
        Console.WriteLine("  dotnet run -- validate-family <familyName>");
        Console.WriteLine("  dotnet run -- validate-all");
        Console.WriteLine("  dotnet run -- zero-forgetting");
        Console.WriteLine("  dotnet run -- trace-mutation <mutationId> [outputPath]");
        Console.WriteLine();
        Console.WriteLine("Example:");
        Console.WriteLine("  dotnet run -- graft P00720 t4_lysozyme PF00959 data/cif/t4_lysozyme_P00720_2LZM.cif");
        Console.WriteLine("  dotnet run -- batch-run P00720 t4-lysozyme 4");
        Console.WriteLine("  dotnet run -- whereami");
        Console.WriteLine("  dotnet run -- trace-mutation MUT_P00720_A_102_LEU");
        Console.WriteLine();
        Console.WriteLine("Environment Variables (Optional):");
        Console.WriteLine("  NEO4J_URI        (default: bolt://localhost:7687)");
        Console.WriteLine("  NEO4J_USERNAME   (default: neo4j)");
        Console.WriteLine("  NEO4J_PASSWORD   (default: mutcert)");
    }

    // =========================================================================
    // Phase 7+8 CLI handlers
    // =========================================================================

    /// <summary>
    /// batch-run: loads all Mutation nodes for a family, feeds them into
    /// MutationAgentPool, and prints a summary when done.
    /// Usage: dotnet run -- batch-run &lt;uniprotId&gt; &lt;familyDbName&gt; [workerCount]
    /// </summary>
    private static async Task<int> RunBatchAsync(string[] args)
    {
        if (args.Length < 3)
        {
            Console.WriteLine("Error: batch-run requires <uniprotId> <familyDbName> [workerCount].");
            return 1;
        }

        string uniprotId    = args[1];
        string familyDbName = args[2];
        int workerCount     = args.Length >= 4 && int.TryParse(args[3], out int wc) ? wc : 4;

        string uri      = Environment.GetEnvironmentVariable("NEO4J_URI")      ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";
        string csvPath  = "data/s2648/s2648.csv";

        Console.WriteLine("================================================================================");
        Console.WriteLine($" MutCert Phase 7 — Batch Run  ({uniprotId} / {familyDbName})");
        Console.WriteLine($" Workers: {workerCount}");
        Console.WriteLine("================================================================================");

        var sw = Stopwatch.StartNew();
        try
        {
            // Load all Mutation nodes from the family database
            var jobs = new System.Collections.Generic.List<MutationJob>();
            using var driver = Neo4j.Driver.GraphDatabase.Driver(uri, AuthTokens.Basic(username, password));
            await using var session = driver.AsyncSession(o => o.WithDatabase(familyDbName));
            await session.ExecuteReadAsync(async tx =>
            {
                var cursor = await tx.RunAsync(@"
                    MATCH (m:Mutation)
                    RETURN m.id AS id, m.mutantType AS mutant, m.seqPos AS pos,
                           m.mutatedChain AS chain");
                while (await cursor.FetchAsync())
                {
                    string mutId   = cursor.Current["id"].As<string>();
                    string mutant  = cursor.Current["mutant"].As<string>();
                    int    seqPos  = cursor.Current["pos"].As<int>();
                    string chain   = cursor.Current["chain"].As<string>();
                    string resId   = $"{uniprotId}_{chain}_{seqPos}";
                    jobs.Add(new MutationJob(uniprotId, familyDbName, resId, mutant, mutId));
                }
            });

            Console.WriteLine($"Loaded {jobs.Count} mutation jobs.");

            await using var pool = new MutationAgentPool(uri, username, password, csvPath, workerCount);
            await pool.StartAsync();

            foreach (var job in jobs)
                await pool.EnqueueAsync(job);

            pool.Complete();
            await pool.DrainAsync();

            sw.Stop();
            int errors    = 0;
            int converged = 0;
            foreach (var r in pool.Results)
            {
                if (r.IsError) errors++;
                else if (r.Converged) converged++;
            }

            Console.WriteLine("================================================================================");
            Console.WriteLine($" Batch complete in {sw.Elapsed.TotalSeconds:F1}s");
            Console.WriteLine($" Total jobs:  {pool.Results.Count}");
            Console.WriteLine($" Converged:   {converged}  ({100.0 * converged / Math.Max(1, pool.Results.Count):F1}%)");
            Console.WriteLine($" Errors:      {errors}");
            Console.WriteLine("================================================================================");
            return 0;
        }
        catch (Exception ex)
        {
            sw.Stop();
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"\nError during batch-run: {ex.Message}");
            Console.WriteLine(ex.StackTrace);
            Console.ResetColor();
            return 2;
        }
    }

    /// <summary>
    /// whereami: prints the FSDE session summary from run-registry.
    /// </summary>
    private static async Task<int> RunWhereAmIAsync()
    {
        string uri      = Environment.GetEnvironmentVariable("NEO4J_URI")      ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

        try
        {
            await using var store = new FsdeSessionStore(uri, username, password);
            string summary = await store.GetSessionSummaryAsync();
            Console.WriteLine(summary);
            return 0;
        }
        catch (Exception ex)
        {
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"Error: {ex.Message}");
            Console.ResetColor();
            return 2;
        }
    }

    /// <summary>
    /// replay: prints the full RunStep trace for a given MutationRun id.
    /// Usage: dotnet run -- replay &lt;runId&gt;
    /// </summary>
    private static async Task<int> RunReplayAsync(string[] args)
    {
        if (args.Length < 2)
        {
            Console.WriteLine("Error: replay requires <runId>.");
            return 1;
        }

        string runId    = args[1];
        string uri      = Environment.GetEnvironmentVariable("NEO4J_URI")      ?? "bolt://localhost:7687";
        string username = Environment.GetEnvironmentVariable("NEO4J_USERNAME") ?? "neo4j";
        string password = Environment.GetEnvironmentVariable("NEO4J_PASSWORD") ?? "mutcert";

        try
        {
            await using var store = new FsdeSessionStore(uri, username, password);
            await store.PrintReplayAsync(runId);
            return 0;
        }
        catch (Exception ex)
        {
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"Error: {ex.Message}");
            Console.ResetColor();
            return 2;
        }
    }
}

