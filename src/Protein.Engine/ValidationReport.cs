using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Protein.Engine;

/// <summary>
/// Generates a self-contained HTML validation report for Phase 9.
///
/// The report includes:
///   - Per-family reliability calibration curve (table + inline bar chart)
///   - Baseline comparison (Mean predictor, MJ-direct, MutCert)
///   - Pass/fail badge for the ≥ 80% @ ±2.0 kcal/mol target
///   - Summary statistics
/// </summary>
public static class ValidationReport
{
    public static void Generate(
        List<FamilyValidationResult> results,
        string outputPath = "data/validation_report.html")
    {
        Directory.CreateDirectory(Path.GetDirectoryName(outputPath)!);
        File.WriteAllText(outputPath, BuildHtml(results), Encoding.UTF8);
        Console.WriteLine($"\n[ValidationReport] Saved → {outputPath}");
    }

    // =========================================================================
    // HTML builder
    // =========================================================================

    private static string BuildHtml(List<FamilyValidationResult> results)
    {
        var sb = new StringBuilder();

        sb.AppendLine("<!DOCTYPE html>");
        sb.AppendLine("<html lang=\"en\">");
        sb.AppendLine(Head());
        sb.AppendLine("<body>");
        sb.AppendLine(Banner(results));

        foreach (var r in results)
        {
            sb.AppendLine(FamilySection(r));
            sb.AppendLine(StructureViewer.Section(r));
        }

        sb.AppendLine(CrossFamilySummary(results));
        sb.AppendLine(Footer());
        sb.AppendLine("</body></html>");
        return sb.ToString();
    }

    // =========================================================================
    // Top-level sections
    // =========================================================================

    private static string Banner(List<FamilyValidationResult> results)
    {
        int totalVal    = results.Sum(r => r.ValCount);
        int totalTrain  = results.Sum(r => r.TrainCount);
        string runAt    = DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss");

        // Overall pass: at least ONE family hits ≥ 80% at ±2.0 kcal/mol
        bool anyPass = results.Any(r =>
            r.CalibrationCurve.FirstOrDefault(c => c.WidthThreshold == 2.0) is { } row &&
            row.MutCertCoverage >= 0.80);

        string badgeCls  = anyPass ? "badge-pass" : "badge-fail";
        string badgeTxt  = anyPass ? "TARGET MET ✓" : "TARGET PENDING";

        return $@"
<div class=""hero"">
  <h1>MutCert — Phase 9 Validation Report</h1>
  <p class=""subtitle"">Reliability Calibration Curve &amp; Baseline Comparison</p>
  <div class=""meta"">
    <span>Run: {runAt}</span>
    <span>Families: {results.Count}</span>
    <span>Train: {totalTrain:N0}  &nbsp; Val: {totalVal:N0}</span>
    <span class=""badge {badgeCls}"">{badgeTxt}</span>
  </div>
  <p class=""target-note"">Phase 1 target: ≥ 80 % coverage at ±2.0 kcal/mol on at least 1 family</p>
</div>";
    }

    private static string FamilySection(FamilyValidationResult r)
    {
        var sb = new StringBuilder();
        sb.AppendLine($"<section class=\"family-section\">");
        sb.AppendLine($"<h2>{r.FamilyName} <span class=\"uniprot\">({r.UniprotId})</span></h2>");

        // Stats row
        sb.AppendLine($@"<div class=""stats-row"">
  <div class=""stat""><span class=""stat-num"">{r.TrainCount:N0}</span><span>training</span></div>
  <div class=""stat""><span class=""stat-num"">{r.ValCount:N0}</span><span>held-out</span></div>
  <div class=""stat""><span class=""stat-num"">{r.Epsilon0:F3}</span><span>ε₀ kcal/mol</span></div>
  <div class=""stat""><span class=""stat-num"">{r.TrainMeanDdg:+0.00;-0.00}</span><span>train mean ΔΔG</span></div>
  <div class=""stat""><span class=""stat-num"">{(double.IsNaN(r.PearsonR) ? "–" : r.PearsonR.ToString("+0.000;-0.000;0.000"))}</span><span>Pearson r</span></div>
  <div class=""stat""><span class=""stat-num"">{(double.IsNaN(r.SpearmanRho) ? "–" : r.SpearmanRho.ToString("+0.000;-0.000;0.000"))}</span><span>Spearman ρ</span></div>
</div>");

        if (r.CalibrationCurve.Count == 0)
        {
            sb.AppendLine("<p class=\"warn\">⚠ Insufficient data to compute calibration curve.</p>");
            sb.AppendLine("</section>");
            return sb.ToString();
        }

        // Convergence rate summary
        var row2 = r.CalibrationCurve.FirstOrDefault(c => c.WidthThreshold == 2.0);
        if (row2 is not null)
        {
            double convRate = r.ValCount > 0 ? (double)row2.MutCertConverged / r.ValCount : 0;
            sb.AppendLine($@"<p class=""conv-rate"">Convergence rate at ≤ 2.0 kcal/mol: 
              <strong>{row2.MutCertConverged}/{r.ValCount} ({convRate*100:F1}%)</strong></p>");
        }

        // Calibration table
        sb.AppendLine("<table class=\"cal-table\">");
        sb.AppendLine(@"<thead><tr>
  <th>Width&nbsp;(kcal/mol)</th>
  <th>MutCert converged</th>
  <th>MutCert coverage</th>
  <th>Recalib coverage</th>
  <th>MJ-direct coverage</th>
  <th>Mean predictor coverage</th>
  <th>Target</th>
</tr></thead><tbody>");;

        foreach (var row in r.CalibrationCurve)
        {
            string targetCell = "";
            if (row.WidthThreshold == 2.0)
            {
                bool pass = row.MutCertCoverage >= 0.80;
                targetCell = pass
                    ? "<span class=\"badge badge-pass\">✓ PASS</span>"
                    : "<span class=\"badge badge-fail\">✗ FAIL</span>";
            }

            sb.AppendLine($@"<tr>
  <td class=""thresh"">±{row.WidthThreshold:F1}</td>
  <td>{row.MutCertConverged}/{row.TotalHeldOut}</td>
  <td>{CoverageBar(row.MutCertCoverage)}</td>  <td>{CoverageBar(row.RecalibCoverage)}</td>  <td>{CoverageBar(row.MjDirectCoverage)}</td>
  <td>{CoverageBar(row.MeanCoverage)}</td>
  <td>{targetCell}</td>
</tr>");
        }
        sb.AppendLine("</tbody></table>");

        // Linear recalibration info
        if (!double.IsNaN(r.LinearSlope))
        {
            string iSign = r.LinearIntercept >= 0 ? "+" : "";
            sb.AppendLine($@"<p class=""note recalib-note"">
  <strong>Linear recalibration (train-fit OLS):</strong>
  corrected = <code>{r.LinearSlope:+0.000;-0.000}</code> × MJ<sub>direct</sub>
  {iSign}<code>{r.LinearIntercept:F3}</code> kcal/mol &nbsp;|
  <em>Recalib cov.</em> re-centres the interval without changing its width.
</p>");
        }

        // Outcome scatter summary
        sb.AppendLine(OutcomeSummary(r));

        sb.AppendLine("</section>");
        return sb.ToString();
    }

    private static string OutcomeSummary(FamilyValidationResult r)
    {
        if (r.Outcomes.Count == 0) return "";

        var converged   = r.Outcomes.Where(o => o.Converged).ToList();
        var notConverged = r.Outcomes.Where(o => !o.Converged).ToList();

        double covCov   = converged.Count > 0
            ? (double)converged.Count(o => o.ExperimentalDdg >= o.FinalLo && o.ExperimentalDdg <= o.FinalHi) / converged.Count
            : 0;

        // Histogram of experimental DDG values (binned)
        var bins = new int[7]; // [-∞,-3],[-3,-2],[-2,-1],[-1,0],[0,1],[1,2],[2,3],[3+]
        double[] edges = [-3, -2, -1, 0, 1, 2, 3];
        foreach (var o in r.Outcomes)
        {
            int b = Array.BinarySearch(edges, o.ExperimentalDdg);
            if (b < 0) b = ~b;
            b = Math.Clamp(b, 0, bins.Length - 1);
            bins[b]++;
        }
        int maxBin = bins.Max();

        var sb = new StringBuilder();
        sb.AppendLine($@"<div class=""outcome-grid"">
<div class=""outcome-card"">
  <div class=""oc-num"">{converged.Count}</div>
  <div class=""oc-lbl"">Converged</div>
</div>
<div class=""outcome-card"">
  <div class=""oc-num"">{notConverged.Count}</div>
  <div class=""oc-lbl"">Not converged</div>
</div>
<div class=""outcome-card"">
  <div class=""oc-num"">{covCov * 100:F1}%</div>
  <div class=""oc-lbl"">Coverage (converged)</div>
</div>
<div class=""outcome-card"">
  <div class=""oc-num"">{(r.Outcomes.Count > 0 ? r.Outcomes.Average(o => Math.Abs(o.ExperimentalDdg - o.FinalDDG)):0):F3}</div>
  <div class=""oc-lbl"">Mean |exp − pred| kcal/mol</div>
</div>
</div>");

        // Mini bar chart for ΔΔG distribution
        string[] labels = ["< −3", "−3 to −2", "−2 to −1", "−1 to 0", "0 to 1", "1 to 2", "> 2"];
        sb.AppendLine("<div class=\"dist-chart\">");
        sb.AppendLine("<h4>Experimental ΔΔG distribution (held-out set)</h4>");
        sb.AppendLine("<div class=\"bars\">");
        for (int i = 0; i < bins.Length && i < labels.Length; i++)
        {
            int pct = maxBin > 0 ? (int)(bins[i] * 100.0 / maxBin) : 0;
            sb.AppendLine($"<div class=\"bar-col\"><div class=\"bar\" style=\"height:{pct}px\" title=\"{bins[i]}\"></div><div class=\"bar-lbl\">{labels[i]}</div></div>");
        }
        sb.AppendLine("</div></div>");

        return sb.ToString();
    }

    private static string CrossFamilySummary(List<FamilyValidationResult> results)
    {
        var sb = new StringBuilder();
        sb.AppendLine("<section class=\"cross-section\">");
        sb.AppendLine("<h2>Cross-family Summary at ±2.0 kcal/mol</h2>");
        sb.AppendLine("<table class=\"cal-table\">");
        sb.AppendLine("<thead><tr><th>Family</th><th>UniProt</th><th>MutCert conv.</th><th>MutCert cov.</th><th>Recalib cov.</th><th>MJ-direct</th><th>Mean pred.</th><th>Pearson r</th><th>Spearman \u03c1</th><th>Pass</th></tr></thead><tbody>");

        foreach (var r in results)
        {
            var row = r.CalibrationCurve.FirstOrDefault(c => c.WidthThreshold == 2.0);
            if (row is null) continue;
            bool pass = row.MutCertCoverage >= 0.80;
            string badge = pass
                ? "<span class=\"badge badge-pass\">✓</span>"
                : "<span class=\"badge badge-fail\">✗</span>";

            string pearsonStr  = double.IsNaN(r.PearsonR)    ? "–" : r.PearsonR.ToString("+0.000;-0.000");
            string spearmanStr = double.IsNaN(r.SpearmanRho) ? "–" : r.SpearmanRho.ToString("+0.000;-0.000");

            sb.AppendLine($@"<tr>
  <td><strong>{r.FamilyName}</strong></td>
  <td>{r.UniprotId}</td>
  <td>{row.MutCertConverged}/{row.TotalHeldOut}</td>
  <td>{CoverageBar(row.MutCertCoverage)}</td>
  <td>{CoverageBar(row.RecalibCoverage)}</td>
  <td>{CoverageBar(row.MjDirectCoverage)}</td>
  <td>{CoverageBar(row.MeanCoverage)}</td>
  <td>{pearsonStr}</td>
  <td>{spearmanStr}</td>
  <td>{badge}</td>
</tr>");
        }
        sb.AppendLine("</tbody></table>");
        sb.AppendLine(@"<p class=""note"">
  Baselines: <strong>MJ-direct</strong> = |ΔΔG<sub>MJ</sub> − exp| ≤ w/2;
  <strong>Mean predictor</strong> = |train_mean − exp| ≤ w/2.
  MutCert coverage denominator is converged mutations only.
</p>");
        sb.AppendLine("</section>");
        return sb.ToString();
    }

    private static string Footer()
        => $"<footer><p>MutCert — Phase 9 Validation Report &nbsp;|&nbsp; {DateTime.Now:yyyy-MM-dd} &nbsp;|&nbsp; Protein.Engine v1.0</p></footer>";

    // =========================================================================
    // Helpers
    // =========================================================================

    private static string CoverageBar(double fraction)
    {
        int pct     = (int)(fraction * 100);
        string color = pct >= 80 ? "#22c55e" : pct >= 60 ? "#f59e0b" : "#ef4444";
        return $@"<div class=""cov-wrap"">
  <div class=""cov-bar"" style=""width:{pct}%;background:{color}""></div>
  <span class=""cov-pct"">{pct}%</span>
</div>";
    }

    // =========================================================================
    // CSS + <head>
    // =========================================================================

    private static string Head() => @"<head>
<meta charset=""UTF-8"">
<meta name=""viewport"" content=""width=device-width, initial-scale=1.0"">
<title>MutCert — Phase 9 Validation Report</title>
<style>
  :root {
    --bg: #0f172a; --surface: #1e293b; --surface2: #273449;
    --border: #334155; --text: #e2e8f0; --muted: #94a3b8;
    --accent: #38bdf8; --pass: #22c55e; --fail: #ef4444;
    --warn: #f59e0b;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { background: var(--bg); color: var(--text); font-family: 'Segoe UI', system-ui, sans-serif; font-size: 14px; }

  /* ── Hero ─────────────────────────────────────────────────────── */
  .hero { background: linear-gradient(135deg, #0c4a6e 0%, #0f172a 100%);
          padding: 2.5rem 2rem 2rem; border-bottom: 1px solid var(--border); }
  .hero h1 { font-size: 1.8rem; color: var(--accent); font-weight: 700; margin-bottom: .35rem; }
  .hero .subtitle { color: var(--muted); margin-bottom: 1rem; }
  .meta { display: flex; gap: 1.5rem; flex-wrap: wrap; margin-bottom: .75rem; }
  .meta span { color: var(--muted); font-size: .85rem; }
  .target-note { color: var(--muted); font-size: .8rem; font-style: italic; }

  /* ── Sections ──────────────────────────────────────────────────── */
  section { padding: 1.75rem 2rem; border-bottom: 1px solid var(--border); }
  h2 { font-size: 1.2rem; margin-bottom: 1rem; color: var(--accent); }
  .uniprot { font-size: .85rem; color: var(--muted); font-weight: 400; }

  /* ── Stats row ─────────────────────────────────────────────────── */
  .stats-row { display: flex; gap: 1rem; flex-wrap: wrap; margin-bottom: 1rem; }
  .stat { background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
          padding: .6rem 1rem; text-align: center; min-width: 90px; }
  .stat-num { display: block; font-size: 1.3rem; font-weight: 700; color: var(--accent); }
  .stat span:last-child { font-size: .75rem; color: var(--muted); }

  /* ── Table ─────────────────────────────────────────────────────── */
  table.cal-table { width: 100%; border-collapse: collapse; margin-bottom: 1.25rem; }
  .cal-table th { background: var(--surface2); color: var(--muted); font-weight: 600;
                  text-align: left; padding: .55rem .75rem; border-bottom: 2px solid var(--border); font-size: .8rem; }
  .cal-table td { padding: .5rem .75rem; border-bottom: 1px solid var(--border); vertical-align: middle; }
  .cal-table tr:hover td { background: var(--surface); }
  .thresh { font-weight: 700; color: var(--accent); }

  /* ── Coverage bar ──────────────────────────────────────────────── */
  .cov-wrap { display: flex; align-items: center; gap: .5rem; min-width: 140px; }
  .cov-bar { height: 10px; border-radius: 5px; min-width: 3px; transition: width .3s; }
  .cov-pct { font-size: .85rem; white-space: nowrap; }

  /* ── Badge ─────────────────────────────────────────────────────── */
  .badge { display: inline-block; padding: .2rem .55rem; border-radius: 4px;
           font-size: .75rem; font-weight: 700; text-transform: uppercase; letter-spacing: .04em; }
  .badge-pass { background: #14532d; color: var(--pass); }
  .badge-fail { background: #450a0a; color: var(--fail); }

  /* ── Outcome grid ──────────────────────────────────────────────── */
  .outcome-grid { display: flex; gap: 1rem; flex-wrap: wrap; margin: 1rem 0; }
  .outcome-card { background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
                  padding: .75rem 1.25rem; text-align: center; min-width: 110px; }
  .oc-num { font-size: 1.4rem; font-weight: 700; color: var(--accent); }
  .oc-lbl { font-size: .72rem; color: var(--muted); margin-top: .2rem; }

  /* ── Mini bar chart ────────────────────────────────────────────── */
  .dist-chart { margin: 1rem 0; }
  .dist-chart h4 { font-size: .8rem; color: var(--muted); margin-bottom: .5rem; }
  .bars { display: flex; align-items: flex-end; gap: 6px; height: 110px; padding-bottom: 24px; position: relative; }
  .bar-col { display: flex; flex-direction: column; align-items: center; flex: 1; }
  .bar { width: 100%; background: var(--accent); border-radius: 3px 3px 0 0; min-height: 2px; }
  .bar-lbl { font-size: .62rem; color: var(--muted); margin-top: 4px; text-align: center; white-space: nowrap; }

  /* ── Cross-family ──────────────────────────────────────────────── */
  .cross-section { background: var(--surface); }
  .conv-rate { font-size: .85rem; color: var(--muted); margin-bottom: .75rem; }
  .note { font-size: .75rem; color: var(--muted); font-style: italic; margin-top: .75rem; }
  .warn { color: var(--warn); padding: .5rem 0; }

  /* ── Footer ────────────────────────────────────────────────────── */
  .recalib-note { border-left: 3px solid var(--warn); padding-left: .75rem; margin-top: .75rem; }

  /* ── 3Dmol viewer ──────────────────────────────────────────────── */
  .structure-section { background: var(--bg); }
  .mol-viewer  { position: relative; width: 100%; height: 460px; margin-top: .5rem;
                 border: 1px solid var(--border); border-radius: 8px; background: #0f172a; }
  .mol-legend  { display: flex; gap: 1rem; flex-wrap: wrap; align-items: center;
                 font-size: .78rem; color: var(--muted); margin: .25rem 0 .5rem; }
  .mol-legend .dot { display: inline-block; width: 10px; height: 10px; border-radius: 50%;
                     margin-right: .35rem; vertical-align: middle; }
  .mol-hint { font-style: italic; }
</style>
<script src=""https://3dmol.org/build/3Dmol-min.js""></script>
</head>";
}
