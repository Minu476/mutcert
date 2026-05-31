using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Channels;
using System.Threading.Tasks;

namespace Protein.Engine;

/// <summary>
/// A single unit of work for the mutation agent pool.
/// </summary>
public sealed record MutationJob(
    string UniprotId,
    string FamilyDbName,
    string ResidueId,
    string MutantType,
    string MutationId);

/// <summary>
/// Outcome reported after a job completes.
/// </summary>
public sealed record MutationJobResult(
    MutationJob Job,
    bool Converged,
    double FinalDDG,
    double FinalLo,
    double FinalHi,
    int TotalSteps,
    Exception? Error = null)
{
    public bool IsError => Error is not null;
}

/// <summary>
/// Channel-based concurrent pool that runs <see cref="MutationAgent"/> workers in parallel.
///
/// Architecture (Phase 7):
/// - A shared <see cref="InMemoryGraph"/> is loaded once per family and stored in a
///   read-only dictionary.  Each worker <b>clones</b> the graph before applying the
///   mutation, so the original is never modified.
/// - Jobs are written to a bounded <see cref="Channel{T}"/> by the producer; N worker
///   Tasks drain it concurrently.
/// - Pause/resume is handled by a <see cref="SemaphoreSlim"/> gate that workers check
///   between jobs.
/// - Certificate writes to <c>run-registry</c> use a separate write transaction
///   (already scoped inside <see cref="ConvergenceSupervisor.ProcessRunAsync"/>).
///
/// Usage:
/// <code>
/// await using var pool = new MutationAgentPool(uri, user, pass, csvPath, workerCount: 4);
/// await pool.StartAsync();
/// await pool.EnqueueAsync(new MutationJob(...));
/// await pool.CompleteAsync();         // signal no more jobs
/// await pool.DrainAsync();            // wait for all workers to finish
/// </code>
/// </summary>
public sealed class MutationAgentPool : IAsyncDisposable
{
    private readonly string _neo4jUri;
    private readonly string _neo4jUser;
    private readonly string _neo4jPass;
    private readonly string _csvPath;
    private readonly int _workerCount;

    private readonly Channel<MutationJob> _channel;
    private readonly SemaphoreSlim _pauseGate = new(1, 1);      // 1 = running, 0 = paused
    private readonly ConcurrentBag<MutationJobResult> _results  = new();

    // Shared read-only graphs, keyed by UniProt ID (loaded lazily on first use).
    private readonly ConcurrentDictionary<string, InMemoryGraph> _sharedGraphs = new();
    private readonly SemaphoreSlim _graphLoadLock = new(1, 1);

    private List<Task>? _workers;
    private CancellationTokenSource? _cts;

    public IReadOnlyCollection<MutationJobResult> Results => _results;

    public MutationAgentPool(
        string neo4jUri,
        string neo4jUser,
        string neo4jPass,
        string csvPath = "data/s2648/s2648.csv",
        int workerCount = 4,
        int channelCapacity = 512)
    {
        if (workerCount < 1)
            throw new ArgumentOutOfRangeException(nameof(workerCount), "Must be ≥ 1.");

        _neo4jUri    = neo4jUri;
        _neo4jUser   = neo4jUser;
        _neo4jPass   = neo4jPass;
        _csvPath     = csvPath;
        _workerCount = workerCount;

        _channel = Channel.CreateBounded<MutationJob>(new BoundedChannelOptions(channelCapacity)
        {
            SingleWriter       = false,
            SingleReader       = false,
            FullMode           = BoundedChannelFullMode.Wait,
            AllowSynchronousContinuations = false
        });
    }

    // -------------------------------------------------------------------------
    // Producer API
    // -------------------------------------------------------------------------

    /// <summary>Enqueues a mutation job. Awaits if the channel is full.</summary>
    public ValueTask EnqueueAsync(MutationJob job, CancellationToken ct = default)
        => _channel.Writer.WriteAsync(job, ct);

    /// <summary>Signals that no more jobs will be enqueued.</summary>
    public void Complete() => _channel.Writer.Complete();

    // -------------------------------------------------------------------------
    // Lifecycle
    // -------------------------------------------------------------------------

    /// <summary>
    /// Starts the worker tasks. Must be called before enqueueing jobs.
    /// </summary>
    public Task StartAsync(CancellationToken ct = default)
    {
        _cts = CancellationTokenSource.CreateLinkedTokenSource(ct);
        _workers = new List<Task>(_workerCount);
        for (int i = 0; i < _workerCount; i++)
        {
            int workerId = i;
            _workers.Add(Task.Run(() => WorkerLoopAsync(workerId, _cts.Token)));
        }
        Console.WriteLine($"[MutationAgentPool] Started {_workerCount} worker(s).");
        return Task.CompletedTask;
    }

    /// <summary>
    /// Waits for all enqueued jobs to be processed.
    /// Call <see cref="Complete"/> before this to signal end of input.
    /// </summary>
    public Task DrainAsync() => _workers is not null
        ? Task.WhenAll(_workers)
        : Task.CompletedTask;

    // -------------------------------------------------------------------------
    // Pause / Resume (FSDE session commands)
    // -------------------------------------------------------------------------

    /// <summary>Pauses the pool between jobs. Currently-running jobs complete first.</summary>
    public async Task PauseAsync()
    {
        if (_pauseGate.CurrentCount > 0)
        {
            await _pauseGate.WaitAsync();
            Console.WriteLine("[MutationAgentPool] Paused.");
        }
    }

    /// <summary>Resumes a paused pool.</summary>
    public void Resume()
    {
        if (_pauseGate.CurrentCount == 0)
        {
            _pauseGate.Release();
            Console.WriteLine("[MutationAgentPool] Resumed.");
        }
    }

    public bool IsPaused => _pauseGate.CurrentCount == 0;

    // -------------------------------------------------------------------------
    // Worker loop
    // -------------------------------------------------------------------------

    private async Task WorkerLoopAsync(int workerId, CancellationToken ct)
    {
        await foreach (var job in _channel.Reader.ReadAllAsync(ct))
        {
            // Check pause gate before starting each job
            await _pauseGate.WaitAsync(ct);
            _pauseGate.Release();

            var result = await ExecuteJobAsync(job, workerId, ct);
            _results.Add(result);

            if (!result.IsError)
            {
                Console.WriteLine(
                    $"[Worker {workerId}] {job.MutationId} | " +
                    $"Converged: {result.Converged} | " +
                    $"DDG: {result.FinalDDG:F3} [{result.FinalLo:F3}, {result.FinalHi:F3}] | " +
                    $"Steps: {result.TotalSteps}");
            }
            else
            {
                Console.WriteLine(
                    $"[Worker {workerId}] ERROR {job.MutationId}: {result.Error!.Message}");
            }
        }
    }

    private async Task<MutationJobResult> ExecuteJobAsync(
        MutationJob job, int workerId, CancellationToken ct)
    {
        try
        {
            // Load shared graph for this family (lazy, thread-safe)
            var sharedGraph = await GetOrLoadGraphAsync(job.UniprotId, ct);

            // Clone: each worker gets its own mutable copy so the shared graph is never
            // modified. This enforces read-only access to the FamilyCluster.
            var graphCopy = sharedGraph.Clone();

            var agent = new MutationAgent(graphCopy);
            var (trace, shellHops) = agent.ApplyMutationAndGetTrace(job.ResidueId, job.MutantType);

            double epsilon0 = EpsilonCalibrator.Calibrate(_csvPath, job.UniprotId);

            // ConvergenceSupervisor uses its own Neo4j session; write transaction is
            // separate from the family-graph read transaction (Phase 7 requirement).
            var supervisor = new ConvergenceSupervisor(
                _neo4jUri, _neo4jUser, _neo4jPass, epsilon0);

            var summary = await supervisor.ProcessRunAsync(
                job.MutationId, job.ResidueId, shellHops, trace);

            return new MutationJobResult(
                job,
                summary.Converged,
                summary.FinalDDG,
                summary.FinalLo,
                summary.FinalHi,
                summary.TotalSteps);
        }
        catch (OperationCanceledException)
        {
            throw;
        }
        catch (Exception ex)
        {
            return new MutationJobResult(
                job,
                Converged: false,
                FinalDDG: 0, FinalLo: 0, FinalHi: 0,
                TotalSteps: 0,
                Error: ex);
        }
    }

    private async Task<InMemoryGraph> GetOrLoadGraphAsync(string uniprotId, CancellationToken ct)
    {
        if (_sharedGraphs.TryGetValue(uniprotId, out var cached))
            return cached;

        await _graphLoadLock.WaitAsync(ct);
        try
        {
            // Double-check after acquiring lock
            if (_sharedGraphs.TryGetValue(uniprotId, out cached))
                return cached;

            Console.WriteLine($"[MutationAgentPool] Loading graph for {uniprotId}...");
            var graph = await GraphLoader.LoadGraphAsync(uniprotId, _neo4jUri, _neo4jUser, _neo4jPass);
            Console.WriteLine($"[MutationAgentPool] Graph loaded: {graph.Nodes.Count} nodes, {graph.Edges.Count} edges.");
            _sharedGraphs[uniprotId] = graph;
            return graph;
        }
        finally
        {
            _graphLoadLock.Release();
        }
    }

    // -------------------------------------------------------------------------
    // IAsyncDisposable
    // -------------------------------------------------------------------------

    public async ValueTask DisposeAsync()
    {
        _cts?.Cancel();
        if (_workers is not null)
        {
            try { await Task.WhenAll(_workers); }
            catch (OperationCanceledException) { }
        }
        _cts?.Dispose();
        _pauseGate.Dispose();
        _graphLoadLock.Dispose();
    }
}
