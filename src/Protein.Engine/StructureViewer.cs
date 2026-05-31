using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace Protein.Engine;

/// <summary>One residue flagged by the wild-type identity audit (Phase 9-C).</summary>
public sealed record ResidueMarker(int SeqPos, string ExpectedWt, string ActualAa, bool Mismatch);

/// <summary>
/// Renders an interactive 3Dmol.js structure viewer section into the Phase 9 HTML report.
///
/// One viewer per family:
///   - cartoon backbone, spectrum-coloured
///   - every validation-site residue marked as a sphere, coloured by
///     |predicted FinalDDG − experimental ΔΔG| (green &lt;1, amber 1–2, red &gt;2)
///   - if audit markers are supplied, wild-type MISMATCH positions are forced to
///     magenta with an "expected X / found Y" hover label
///
/// Zero remote structure fetch: the local mmCIF the engine already parsed is embedded
/// inline as a text/plain block, so CI2/barnase (AlphaFold models with no PDB id)
/// render the same way T4 lysozyme (2LZM) does — from the exact files the pipeline read.
///
/// ─────────────────────────────────────────────────────────────────────────────
/// INTEGRATION (deferred — do after CI2/barnase sequence offset fix, Phase 9-D)
/// ─────────────────────────────────────────────────────────────────────────────
///
/// 1. ValidationReport.cs — BuildHtml loop, add after FamilySection():
///
///      foreach (var r in results)
///      {
///          sb.AppendLine(FamilySection(r));
///          sb.AppendLine(StructureViewer.Section(r));          // basic
///          // or, once AuditMarkers is wired:
///          // sb.AppendLine(StructureViewer.Section(r, r.AuditMarkers));
///      }
///
/// 2. ValidationReport.cs — Head(), inside &lt;style&gt; before closing &lt;/style&gt;:
///
///      /* ── 3Dmol viewer ──────────────────────────────────────────── */
///      .structure-section { background: var(--bg); }
///      .mol-viewer  { position: relative; width: 100%; height: 460px; margin-top: .5rem;
///                     border: 1px solid var(--border); border-radius: 8px; background: #0f172a; }
///      .mol-legend  { display: flex; gap: 1rem; flex-wrap: wrap; align-items: center;
///                     font-size: .78rem; color: var(--muted); margin: .25rem 0 .5rem; }
///      .mol-legend .dot { display: inline-block; width: 10px; height: 10px; border-radius: 50%;
///                         margin-right: .35rem; vertical-align: middle; }
///      .mol-hint { font-style: italic; }
///
/// 3. ValidationReport.cs — Head(), just before &lt;/head&gt;:
///
///      &lt;script src="https://3dmol.org/build/3Dmol-min.js"&gt;&lt;/script&gt;
///      (or pin a version: https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.0/3Dmol-min.js)
///
/// 4. Optional — thread AuditMarkers through FamilyValidationResult:
///    a. Add trailing optional field to FamilyValidationResult record:
///          IReadOnlyList&lt;ResidueMarker&gt;? AuditMarkers = null);
///    b. Build auditMarkers in RunFamilyAsync (after WT-audit step, step 3.5):
///          var auditMarkers = val
///              .GroupBy(m => m.ResidueId).Select(g => g.First())
///              .Select(m => {
///                  int pos = int.TryParse(m.ResidueId.Split('_').Last(), out var p) ? p : 0;
///                  string actual = frozenGraph.Nodes.TryGetValue(m.ResidueId, out var n) ? n.ResidueName : "—";
///                  bool mismatch = !string.Equals(actual, m.WildType3, StringComparison.OrdinalIgnoreCase);
///                  return new ResidueMarker(pos, m.WildType3, actual, mismatch);
///              }).ToList();
///    c. Pass auditMarkers as last argument to new FamilyValidationResult(...)
///    d. Switch call to: StructureViewer.Section(r, r.AuditMarkers)
///
/// NOTE: Until the CI2/barnase sequence offset bug is fixed (Phase 9-D), spheres
/// on those families will sit on the wrong residues and magenta flags will appear
/// across the structure. That is expected — it makes the misalignment visible.
/// After the offset fix, re-run and watch the spheres relocate and magenta clear.
/// </summary>
public static class StructureViewer
{
    private const string CifDir = "data/cif";

    public static string Section(FamilyValidationResult r, IReadOnlyList<ResidueMarker>? markers = null)
    {
        string? cifPath = ResolveCif(r.UniprotId);
        if (cifPath is null)
            return $"<section class=\"family-section structure-section\"><h2>3D structure — {r.FamilyName}</h2>" +
                   $"<p class=\"warn\">⚠ No mmCIF found in {CifDir} matching {r.UniprotId} — 3D viewer skipped.</p></section>";

        string cifRaw;
        try { cifRaw = File.ReadAllText(cifPath); }
        catch (Exception ex)
        {
            return $"<section class=\"family-section structure-section\"><h2>3D structure — {r.FamilyName}</h2>" +
                   $"<p class=\"warn\">⚠ Could not read {cifPath}: {ex.Message}</p></section>";
        }

        var sites = AggregateSites(r);
        var markerByPos = (markers ?? new List<ResidueMarker>())
            .GroupBy(m => m.SeqPos)
            .ToDictionary(g => g.Key, g => g.First());
        int mismatchCount = markerByPos.Values.Count(m => m.Mismatch);

        string viewerId = "mol_" + Safe(r.UniprotId);
        string cifId    = "cif_" + Safe(r.UniprotId);

        // Build the JS site array literal: {resi, chain, color, tip}
        var jsSites = new StringBuilder();
        foreach (var s in sites)
        {
            bool mism = markerByPos.TryGetValue(s.SeqPos, out var mk) && mk.Mismatch;
            string color = mism ? "#d946ef" : ErrColor(s.MeanAbsErr);
            string tip = mism
                ? $"resi {mk!.SeqPos}: expected {mk.ExpectedWt}, found {mk.ActualAa} — WT MISMATCH"
                : $"resi {s.SeqPos}: |err| {s.MeanAbsErr:F2} kcal/mol (n={s.Count})";

            jsSites.Append("{resi:").Append(s.SeqPos)
                   .Append(",chain:").Append(JsStr(s.Chain))
                   .Append(",color:'").Append(color).Append('\'')
                   .Append(",tip:").Append(JsStr(tip)).Append("},");
        }

        string js = JsTemplate
            .Replace("%VIEWER_ID%", viewerId)
            .Replace("%CIF_ID%", cifId)
            .Replace("%SITES%", jsSites.ToString());

        string mmNote = markers is null
            ? "<span class=\"mol-hint\">(pass audit markers to flag WT mismatches in magenta)</span>"
            : (mismatchCount > 0
                ? $"<span style=\"color:#d946ef;font-weight:700\">{mismatchCount} wild-type mismatch position(s) flagged</span>"
                : "<span style=\"color:#22c55e\">all audited positions match the structure</span>");

        var sb = new StringBuilder();
        sb.AppendLine("<section class=\"family-section structure-section\">");
        sb.AppendLine($"<h2>3D structure — {r.FamilyName} <span class=\"uniprot\">({Path.GetFileName(cifPath)})</span></h2>");
        sb.AppendLine($@"<div class=""mol-legend"">
  <span><i class=""dot"" style=""background:#22c55e""></i>|err| &lt; 1</span>
  <span><i class=""dot"" style=""background:#f59e0b""></i>1–2</span>
  <span><i class=""dot"" style=""background:#ef4444""></i>&gt; 2 kcal/mol</span>
  <span><i class=""dot"" style=""background:#d946ef""></i>WT mismatch</span>
  {mmNote}
  <span class=""mol-hint"">drag to rotate · scroll to zoom · hover a sphere for details</span>
</div>");
        sb.AppendLine($"<div id=\"{viewerId}\" class=\"mol-viewer\"></div>");

        // Embed mmCIF as text/plain so it neither executes nor needs JS escaping.
        sb.AppendLine($"<script type=\"text/plain\" id=\"{cifId}\">");
        sb.Append(cifRaw);
        sb.AppendLine();
        sb.AppendLine("</script>");

        sb.AppendLine("<script>");
        sb.AppendLine(js);
        sb.AppendLine("</script>");
        sb.AppendLine("</section>");
        return sb.ToString();
    }

    // ------------------------------------------------------------------ helpers

    private static string? ResolveCif(string uniprotId)
    {
        if (!Directory.Exists(CifDir)) return null;
        return Directory.GetFiles(CifDir, "*.cif")
            .FirstOrDefault(f => Path.GetFileName(f).Contains(uniprotId, StringComparison.OrdinalIgnoreCase));
    }

    private readonly record struct SiteAgg(int SeqPos, string Chain, double MeanAbsErr, int Count);

    private static List<SiteAgg> AggregateSites(FamilyValidationResult r)
    {
        // Multiple mutations can hit the same position; aggregate to mean |error|.
        var acc = new Dictionary<int, (string chain, double sum, int n)>();
        foreach (var o in r.Outcomes)
        {
            if (double.IsNaN(o.ExperimentalDdg)) continue;
            var (chain, pos) = ParseChainPos(o.MutationId);
            if (pos is null) continue;
            double err = Math.Abs(o.FinalDDG - o.ExperimentalDdg);
            if (acc.TryGetValue(pos.Value, out var cur))
                acc[pos.Value] = (cur.chain, cur.sum + err, cur.n + 1);
            else
                acc[pos.Value] = (chain, err, 1);
        }
        return acc
            .Select(kv => new SiteAgg(kv.Key, kv.Value.chain, kv.Value.sum / kv.Value.n, kv.Value.n))
            .OrderBy(s => s.SeqPos)
            .ToList();
    }

    // MutationId format: MUT_{uniprot}_{chain}_{seqPos}_{mutant3}
    private static (string chain, int? pos) ParseChainPos(string mutationId)
    {
        var p = mutationId.Split('_');
        if (p.Length < 5) return ("", null);
        return int.TryParse(p[3], NumberStyles.Integer, CultureInfo.InvariantCulture, out int pos)
            ? (p[2], pos)
            : ("", null);
    }

    private static string ErrColor(double e) => e < 1.0 ? "#22c55e" : e <= 2.0 ? "#f59e0b" : "#ef4444";

    private static string Safe(string s) => new string(s.Select(c => char.IsLetterOrDigit(c) ? c : '_').ToArray());

    private static string JsStr(string s)
    {
        var b = new StringBuilder("'");
        foreach (char c in s)
        {
            if (c is '\\' or '\'') b.Append('\\').Append(c);
            else if (c is '\n' or '\r') b.Append(' ');
            else b.Append(c);
        }
        return b.Append('\'').ToString();
    }

    // Token-substituted (not interpolated) so literal { } braces survive.
    private const string JsTemplate = @"
(function () {
  function draw() {
    if (typeof $3Dmol === 'undefined') { setTimeout(draw, 120); return; }   // wait for CDN
    var el = document.getElementById('%VIEWER_ID%');
    if (!el) return;
    var cif = document.getElementById('%CIF_ID%').textContent;
    var viewer = $3Dmol.createViewer(el, { backgroundColor: '#0f172a' });
    viewer.addModel(cif, 'cif');
    viewer.setStyle({}, { cartoon: { color: 'spectrum' } });

    var sites = [%SITES%];
    sites.forEach(function (s) {
      var sel = { resi: s.resi };
      if (s.chain) sel.chain = s.chain;
      viewer.addStyle(sel, { stick: { radius: 0.15 } });
      viewer.addStyle(sel, { sphere: { scale: 0.35, color: s.color } });
    });

    viewer.setHoverable({}, true,
      function (atom, vw) {
        if (!atom || atom.__lbl) return;
        var hit = null;
        for (var i = 0; i < sites.length; i++) {
          if (sites[i].resi === atom.resi && (!sites[i].chain || sites[i].chain === atom.chain)) { hit = sites[i]; break; }
        }
        var txt = hit ? hit.tip : (atom.resn + ' ' + atom.resi);
        atom.__lbl = vw.addLabel(txt, { position: atom, backgroundColor: '#1e293b', fontColor: '#e2e8f0', fontSize: 12, borderColor: '#334155', borderThickness: 1 });
        vw.render();
      },
      function (atom, vw) { if (atom && atom.__lbl) { vw.removeLabel(atom.__lbl); delete atom.__lbl; vw.render(); } });

    viewer.zoomTo();
    viewer.render();
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', draw); else draw();
})();
";
}
