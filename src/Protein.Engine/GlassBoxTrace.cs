using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Neo4j.Driver;

namespace Protein.Engine;

/// <summary>
/// Phase 11 — Glass-box trace generator.
///
/// Queries the run-registry for a stored MutationRun and its RunStep chain,
/// then renders a self-contained HTML showing the full causal chain from
/// mutation to convergence certificate — every BFS step visible and annotated.
///
/// CLI: dotnet run -- trace-mutation &lt;mutationId&gt; [outputPath]
/// Output: data/trace_{mutationId}.html (default)
/// </summary>
public sealed class GlassBoxTrace : IAsyncDisposable
{
    private readonly IDriver _driver;

    public GlassBoxTrace(string uri, string username, string password)
    {
        _driver = GraphDatabase.Driver(uri, AuthTokens.Basic(username, password));
    }

    // =========================================================================
    // Public entry point
    // =========================================================================

    public async Task<string> GenerateAsync(string mutationId, string outputPath)
    {
        var run   = await LoadRunAsync(mutationId);
        var steps = await LoadStepsAsync(run.RunId);

        if (steps.Count == 0)
            throw new InvalidOperationException(
                $"No RunStep records found for mutationId '{mutationId}'. " +
                "Run validate-all first to populate run-registry.");

        string html = BuildHtml(run, steps);
        Directory.CreateDirectory(Path.GetDirectoryName(outputPath)!);
        File.WriteAllText(outputPath, html, Encoding.UTF8);
        return outputPath;
    }

    // =========================================================================
    // Neo4j queries
    // =========================================================================

    private sealed record RunInfo(
        string RunId, string MutationId, string TargetResidueId,
        double FinalDDG, double FinalLo, double FinalHi,
        bool Converged, int TotalSteps, double Epsilon0,
        double ExperimentalDdg, long Timestamp);

    private async Task<RunInfo> LoadRunAsync(string mutationId)
    {
        await using var session = _driver.AsyncSession(o => o.WithDatabase("run-registry"));
        return await session.ExecuteReadAsync(async tx =>
        {
            var cursor = await tx.RunAsync(@"
                MATCH (m:MutationRun { mutationId: $mutationId })
                RETURN m.id AS runId, m.mutationId AS mutId,
                       m.targetResidueId AS resId,
                       m.finalDDG AS ddg, m.finalLo AS lo, m.finalHi AS hi,
                       m.converged AS conv, m.totalSteps AS steps,
                       m.epsilon0 AS eps, m.experimentalDdg AS expDdg,
                       m.timestamp AS ts
                ORDER BY m.timestamp DESC
                LIMIT 1",
                new { mutationId });

            if (!await cursor.FetchAsync())
                throw new InvalidOperationException(
                    $"MutationRun with mutationId '{mutationId}' not found in run-registry.");

            var r = cursor.Current;
            return new RunInfo(
                r["runId"].As<string>(),
                r["mutId"].As<string>(),
                r["resId"].As<string>(),
                r["ddg"].As<double>(),
                r["lo"].As<double>(),
                r["hi"].As<double>(),
                r["conv"].As<bool>(),
                r["steps"].As<int>(),
                r["eps"].As<double>(),
                r["expDdg"] is null || r["expDdg"] == null
                    ? double.NaN
                    : r["expDdg"].As<double>(),
                r["ts"].As<long>());
        });
    }

    private async Task<List<RunStepRecord>> LoadStepsAsync(string runId)
    {
        await using var session = _driver.AsyncSession(o => o.WithDatabase("run-registry"));
        return await session.ExecuteReadAsync(async tx =>
        {
            var result = new List<RunStepRecord>();
            var cursor = await tx.RunAsync(@"
                MATCH (m:MutationRun { id: $runId })-[:HAS_STEP]->(s:RunStep)
                RETURN s.stepIndex AS idx, s.nodeId AS node, s.deltaE AS delta,
                       s.hopDist AS hop, s.lo AS lo, s.hi AS hi
                ORDER BY s.stepIndex",
                new { runId });
            while (await cursor.FetchAsync())
                result.Add(new RunStepRecord(
                    cursor.Current["idx"].As<int>(),
                    cursor.Current["node"].As<string>(),
                    cursor.Current["delta"].As<double>(),
                    cursor.Current["hop"].As<int>(),
                    cursor.Current["lo"].As<double>(),
                    cursor.Current["hi"].As<double>()));
            return result;
        });
    }

    // =========================================================================
    // HTML rendering
    // =========================================================================

    private static string BuildHtml(RunInfo run, List<RunStepRecord> steps)
    {
        var sb = new StringBuilder();
        sb.AppendLine("<!DOCTYPE html>");
        sb.AppendLine("<html lang=\"en\">");
        sb.AppendLine(Head(run));
        sb.AppendLine("<body>");
        sb.AppendLine(Hero(run));
        sb.AppendLine(CertificateBox(run));
        sb.AppendLine(NarrowingChart(run, steps));
        sb.AppendLine(StepTable(steps));
        sb.AppendLine(Footer());
        sb.AppendLine("</body></html>");
        return sb.ToString();
    }

    private static string Head(RunInfo run) => $@"<head>
<meta charset=""UTF-8"">
<meta name=""viewport"" content=""width=device-width, initial-scale=1.0"">
<title>MutCert Glass-Box Trace — {run.MutationId}</title>
<style>
  :root {{
    --bg: #0f172a; --surface: #1e293b; --surface2: #273449;
    --border: #334155; --text: #e2e8f0; --muted: #94a3b8;
    --accent: #38bdf8; --pass: #22c55e; --fail: #ef4444;
    --warn: #f59e0b; --hop0: #818cf8; --hop1: #38bdf8; --hop2: #34d399; --hop3: #fb923c;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: var(--bg); color: var(--text);
          font-family: 'Segoe UI', system-ui, sans-serif; font-size: 14px; }}

  /* Hero */
  .hero {{ background: linear-gradient(135deg, #0c4a6e 0%, #0f172a 100%);
           padding: 2.5rem 2rem 2rem; border-bottom: 1px solid var(--border); }}
  .hero h1 {{ font-size: 1.6rem; color: var(--accent); font-weight: 700; margin-bottom: .3rem; }}
  .hero .sub {{ color: var(--muted); font-size: .9rem; margin-bottom: 1rem; }}
  .kv-grid {{ display: flex; gap: 1.5rem; flex-wrap: wrap; }}
  .kv {{ background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
         padding: .55rem 1rem; text-align: center; min-width: 110px; }}
  .kv-val {{ display: block; font-size: 1.2rem; font-weight: 700; color: var(--accent); }}
  .kv-lbl {{ font-size: .72rem; color: var(--muted); }}

  /* Certificate box */
  .cert-box {{ margin: 1.5rem 2rem; padding: 1.25rem 1.5rem;
               border-radius: 10px; border: 2px solid; }}
  .cert-box.issued {{ border-color: var(--pass); background: rgba(34,197,94,.07); }}
  .cert-box.not-issued {{ border-color: var(--fail); background: rgba(239,68,68,.07); }}
  .cert-box h2 {{ font-size: 1rem; margin-bottom: .75rem; }}
  .cert-box.issued h2 {{ color: var(--pass); }}
  .cert-box.not-issued h2 {{ color: var(--fail); }}
  .cert-conditions {{ list-style: none; display: flex; flex-direction: column; gap: .35rem; }}
  .cert-conditions li {{ font-size: .85rem; }}
  .cert-conditions li::before {{ margin-right: .5rem; }}
  .cond-pass::before {{ content: '✓'; color: var(--pass); }}
  .cond-fail::before {{ content: '✗'; color: var(--fail); }}

  /* Chart */
  .chart-section {{ margin: 0 2rem 1.5rem; }}
  .chart-section h2 {{ font-size: 1rem; color: var(--accent); margin-bottom: .75rem; }}
  svg.narrowing {{ width: 100%; height: 220px; background: var(--surface);
                   border: 1px solid var(--border); border-radius: 8px; }}

  /* Step table */
  .table-section {{ padding: 0 2rem 2rem; }}
  .table-section h2 {{ font-size: 1rem; color: var(--accent); margin-bottom: .75rem; }}
  table.steps {{ width: 100%; border-collapse: collapse; font-size: .8rem; }}
  .steps th {{ background: var(--surface2); color: var(--muted); font-weight: 600;
               text-align: left; padding: .4rem .65rem; border-bottom: 2px solid var(--border);
               position: sticky; top: 0; }}
  .steps td {{ padding: .35rem .65rem; border-bottom: 1px solid var(--border); }}
  .steps tr:hover td {{ background: var(--surface); }}
  .hop-0 {{ color: var(--hop0); font-weight: 700; }}
  .hop-1 {{ color: var(--hop1); }}
  .hop-2 {{ color: var(--hop2); }}
  .hop-3 {{ color: var(--hop3); }}
  .neg {{ color: #f87171; }}
  .pos {{ color: #4ade80; }}
  .interval {{ font-family: 'Courier New', monospace; font-size: .78rem; color: var(--muted); }}
  .width-narrow {{ color: var(--pass); }}
  .width-wide {{ color: var(--fail); }}

  /* Footer */
  footer {{ padding: 1rem 2rem; color: var(--muted); font-size: .75rem;
            border-top: 1px solid var(--border); text-align: center; }}
</style>
</head>";

    private static string Hero(RunInfo run)
    {
        bool expKnown = !double.IsNaN(run.ExperimentalDdg);
        string expStr = expKnown ? $"{run.ExperimentalDdg:+0.000;-0.000} kcal/mol" : "—";
        string covered = expKnown
            ? (run.ExperimentalDdg >= run.FinalLo && run.ExperimentalDdg <= run.FinalHi
               ? "✓ inside" : "✗ outside")
            : "—";
        string runAt = DateTimeOffset.FromUnixTimeMilliseconds(run.Timestamp).LocalDateTime
                                     .ToString("yyyy-MM-dd HH:mm:ss");

        string[] parts = run.TargetResidueId.Split('_');
        string family = parts.Length >= 2 ? parts[0] : run.TargetResidueId;

        return $@"<div class=""hero"">
  <h1>Glass-Box Causal Trace</h1>
  <p class=""sub"">Mutation: <strong>{run.MutationId}</strong> &nbsp;|&nbsp; Family: {family} &nbsp;|&nbsp; Run: {runAt}</p>
  <div class=""kv-grid"">
    <div class=""kv""><span class=""kv-val"">{run.TotalSteps}</span><span class=""kv-lbl"">BFS steps</span></div>
    <div class=""kv""><span class=""kv-val"">{run.FinalDDG:+0.000;-0.000}</span><span class=""kv-lbl"">predicted ΔΔG (kcal/mol)</span></div>
    <div class=""kv""><span class=""kv-val"">{run.FinalHi - run.FinalLo:F3}</span><span class=""kv-lbl"">final width (kcal/mol)</span></div>
    <div class=""kv""><span class=""kv-val"">[{run.FinalLo:F2}, {run.FinalHi:F2}]</span><span class=""kv-lbl"">final interval</span></div>
    <div class=""kv""><span class=""kv-val"">{run.Epsilon0:F4}</span><span class=""kv-lbl"">ε₀ (kcal/mol)</span></div>
    <div class=""kv""><span class=""kv-val"">{expStr}</span><span class=""kv-lbl"">experimental ΔΔG</span></div>
    <div class=""kv""><span class=""kv-val"">{covered}</span><span class=""kv-lbl"">exp. inside interval</span></div>
  </div>
</div>";
    }

    private static string CertificateBox(RunInfo run)
    {
        double width   = run.FinalHi - run.FinalLo;
        bool c1 = width <= 2.0;
        bool c2 = run.Converged;    // monotone tracked by supervisor
        bool c3 = run.TotalSteps >= 5;
        bool issued = run.Converged;

        string cls  = issued ? "issued" : "not-issued";
        string head = issued ? "⬛ Convergence Certificate — ISSUED" : "⬛ Convergence Certificate — NOT ISSUED";

        string Cond(bool pass, string text) =>
            $"<li class=\"cond-{(pass ? "pass" : "fail")}\">{text}</li>";

        return $@"<div class=""cert-box {cls}"">
  <h2>{head}</h2>
  <p style=""font-size:.8rem;color:var(--muted);margin-bottom:.75rem"">
    Formal conditions (spec §6, Phase 1 threshold):
  </p>
  <ul class=""cert-conditions"">
    {Cond(c1, $"C1: width(k) = {width:F3} kcal/mol ≤ 2.0 kcal/mol (Phase 1 threshold)")}
    {Cond(c2, "C2: interval width non-increasing for ≥ 5 consecutive steps (monotone narrowing)")}
    {Cond(c3, $"C3: ≥ 5 BFS steps completed (8 Å shell evaluation) — {run.TotalSteps} steps")}
  </ul>
  <p style=""font-size:.75rem;color:var(--muted);margin-top:.75rem;font-style:italic"">
    Phase 2 target: tighten threshold to 1.0 kcal/mol and extend shell to 12 Å.
  </p>
</div>";
    }

    private static string NarrowingChart(RunInfo run, List<RunStepRecord> steps)
    {
        // SVG inline chart: X = step index, Y = [lo(k), hi(k)] interval band
        int n      = steps.Count;
        if (n == 0) return "";

        double allLo  = steps.Min(s => s.Lo);
        double allHi  = steps.Max(s => s.Hi);
        double range  = allHi - allLo;
        if (range <= 0) range = 1.0;

        const int W = 900, H = 180, PAD_L = 55, PAD_R = 20, PAD_T = 15, PAD_B = 25;
        int plotW = W - PAD_L - PAD_R;
        int plotH = H - PAD_T - PAD_B;

        double Px(int i) => PAD_L + (double)i / Math.Max(n - 1, 1) * plotW;
        double Py(double v) => PAD_T + plotH - (v - allLo) / range * plotH;

        // Build the shaded band polygon (lo path top, hi path reversed)
        var loPoints = steps.Select((s, i) => $"{Px(i).ToString("F1", CultureInfo.InvariantCulture)},{Py(s.Lo).ToString("F1", CultureInfo.InvariantCulture)}");
        var hiPoints = steps.Select((s, i) => $"{Px(i).ToString("F1", CultureInfo.InvariantCulture)},{Py(s.Hi).ToString("F1", CultureInfo.InvariantCulture)}").Reverse();
        string band = string.Join(" ", loPoints.Concat(hiPoints));

        // Mid-line (running DDG)
        string midLine = string.Join(" ", steps.Select((s, i) =>
            $"{(i == 0 ? "M" : "L")}{Px(i).ToString("F1", CultureInfo.InvariantCulture)},{Py((s.Lo + s.Hi) / 2.0).ToString("F1", CultureInfo.InvariantCulture)}"));

        // Experimental DDG line (if known)
        string expLine = "";
        if (!double.IsNaN(run.ExperimentalDdg) && run.ExperimentalDdg >= allLo && run.ExperimentalDdg <= allHi)
        {
            double ey = Py(run.ExperimentalDdg);
            expLine = $@"<line x1=""{PAD_L}"" y1=""{ey.ToString("F1", CultureInfo.InvariantCulture)}"" x2=""{W - PAD_R}"" y2=""{ey.ToString("F1", CultureInfo.InvariantCulture)}"" stroke=""#f59e0b"" stroke-width=""1.5"" stroke-dasharray=""6,3"" opacity="".85""/>";
        }

        // Y-axis ticks
        var sb = new StringBuilder();
        int nTicks = 4;
        for (int t = 0; t <= nTicks; t++)
        {
            double v = allLo + range * t / nTicks;
            double y = Py(v);
            sb.AppendLine($"<line x1=\"{PAD_L}\" y1=\"{y.ToString("F1", CultureInfo.InvariantCulture)}\" x2=\"{W - PAD_R}\" y2=\"{y.ToString("F1", CultureInfo.InvariantCulture)}\" stroke=\"#334155\" stroke-width=\"1\"/>");
            sb.AppendLine($"<text x=\"{PAD_L - 4}\" y=\"{(y + 4).ToString("F1", CultureInfo.InvariantCulture)}\" text-anchor=\"end\" font-size=\"9\" fill=\"#94a3b8\">{v:F1}</text>");
        }

        return $@"<div class=""chart-section"">
  <h2>Interval Narrowing — [lo(k), hi(k)] over BFS steps</h2>
  <svg class=""narrowing"" viewBox=""0 0 {W} {H}"" preserveAspectRatio=""xMidYMid meet"">
    <!-- Grid -->
    {sb}
    <!-- Shaded interval band -->
    <polygon points=""{band}"" fill=""#38bdf8"" opacity="".15""/>
    <!-- lo and hi boundary lines -->
    <polyline points=""{string.Join(" ", steps.Select((s, i) => $"{Px(i).ToString("F1", CultureInfo.InvariantCulture)},{Py(s.Lo).ToString("F1", CultureInfo.InvariantCulture)}"))}"" fill=""none"" stroke=""#38bdf8"" stroke-width=""1.5"" opacity="".7""/>
    <polyline points=""{string.Join(" ", steps.Select((s, i) => $"{Px(i).ToString("F1", CultureInfo.InvariantCulture)},{Py(s.Hi).ToString("F1", CultureInfo.InvariantCulture)}"))}"" fill=""none"" stroke=""#38bdf8"" stroke-width=""1.5"" opacity="".7""/>
    <!-- Running DDG midpoint -->
    <path d=""{midLine}"" fill=""none"" stroke=""#e2e8f0"" stroke-width=""1.5""/>
    <!-- Experimental DDG -->
    {expLine}
    <!-- Axes -->
    <line x1=""{PAD_L}"" y1=""{PAD_T}"" x2=""{PAD_L}"" y2=""{H - PAD_B}"" stroke=""#475569"" stroke-width=""1""/>
    <line x1=""{PAD_L}"" y1=""{H - PAD_B}"" x2=""{W - PAD_R}"" y2=""{H - PAD_B}"" stroke=""#475569"" stroke-width=""1""/>
    <text x=""{W / 2}"" y=""{H - 5}"" text-anchor=""middle"" font-size=""10"" fill=""#94a3b8"">BFS step</text>
    <text x=""10"" y=""{(H / 2)}"" text-anchor=""middle"" font-size=""10"" fill=""#94a3b8"" transform=""rotate(-90 10 {H / 2})"">ΔΔG (kcal/mol)</text>
    <!-- Legend -->
    <rect x=""{W - PAD_R - 160}"" y=""{PAD_T}"" width=""145"" height=""58"" rx=""4"" fill=""#1e293b"" stroke=""#334155""/>
    <rect x=""{W - PAD_R - 150}"" y=""{PAD_T + 8}"" width=""14"" height=""8"" fill=""#38bdf8"" opacity="".35""/>
    <text x=""{W - PAD_R - 132}"" y=""{PAD_T + 16}"" font-size=""9"" fill=""#94a3b8"">Interval [lo, hi]</text>
    <line x1=""{W - PAD_R - 150}"" y1=""{PAD_T + 26}"" x2=""{W - PAD_R - 136}"" y2=""{PAD_T + 26}"" stroke=""#e2e8f0"" stroke-width=""1.5""/>
    <text x=""{W - PAD_R - 132}"" y=""{PAD_T + 30}"" font-size=""9"" fill=""#94a3b8"">Running ΔΔG</text>
    <line x1=""{W - PAD_R - 150}"" y1=""{PAD_T + 44}"" x2=""{W - PAD_R - 136}"" y2=""{PAD_T + 44}"" stroke=""#f59e0b"" stroke-width=""1.5"" stroke-dasharray=""4,2""/>
    <text x=""{W - PAD_R - 132}"" y=""{PAD_T + 48}"" font-size=""9"" fill=""#94a3b8"">Experimental ΔΔG</text>
  </svg>
</div>";
    }

    private static string StepTable(List<RunStepRecord> steps)
    {
        var sb = new StringBuilder();
        sb.AppendLine("<div class=\"table-section\">");
        sb.AppendLine("<h2>BFS Step-by-Step Causal Chain</h2>");
        sb.AppendLine("<table class=\"steps\">");
        sb.AppendLine(@"<thead><tr>
  <th>Step</th>
  <th>Residue node</th>
  <th>Hop</th>
  <th>ΔE (kcal/mol)</th>
  <th>Cumulative ΔΔG</th>
  <th>Interval [lo, hi]</th>
  <th>Width</th>
</tr></thead><tbody>");

        double cumDDG = 0;
        foreach (var s in steps)
        {
            cumDDG += s.DeltaEKcal;
            string hopCls = s.HopDistance switch
            {
                0 => "hop-0", 1 => "hop-1", 2 => "hop-2", _ => "hop-3"
            };
            string deltaCls = s.DeltaEKcal >= 0 ? "pos" : "neg";
            double w = s.Hi - s.Lo;
            string wCls = w <= 2.0 ? "width-narrow" : "width-wide";

            sb.AppendLine($@"<tr>
  <td>{s.StepIndex}</td>
  <td><code>{s.NodeId}</code></td>
  <td class=""{hopCls}"">+{s.HopDistance}</td>
  <td class=""{deltaCls}"">{s.DeltaEKcal:+0.000;-0.000}</td>
  <td>{cumDDG:+0.000;-0.000}</td>
  <td class=""interval"">[{s.Lo:F3}, {s.Hi:F3}]</td>
  <td class=""{wCls}"">{w:F3}</td>
</tr>");
        }

        sb.AppendLine("</tbody></table>");
        sb.AppendLine("</div>");
        return sb.ToString();
    }

    private static string Footer() => $@"
<footer>
  MutCert — Glass-Box Trace &nbsp;|&nbsp;
  Generated {DateTime.Now:yyyy-MM-dd HH:mm:ss} &nbsp;|&nbsp;
  Phase 11 — Patent &amp; Paper Prep
</footer>";

    public async ValueTask DisposeAsync()
    {
        await _driver.DisposeAsync();
    }
}
