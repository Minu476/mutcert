# MutCert — Differentiation vs State-of-the-Art

**Phase 11 Preprint Prep | MutCert v1.0**  
**Keywords:** differentiation, mutpred, esm, foldx, rosetta, ddgun, glass-box, transparency, ratchet-heuristic, phase11

---

## 1. Positioning Statement

MutCert adds a **transparency layer** on top of a statistical contact-energy calculation. The two concrete contributions are:

1. **Glass-box BFS trace**: every residue visited, its hop distance, its energy contribution, and the resulting interval — all stored and inspectable. No existing method provides this.
2. **Convergence-tracking heuristic**: a running interval that narrows monotonically via the max/min ratchet and emits a record when precision conditions C1–C3 are met.

The claim is *not* "certified measurement protocol" — that framing puts the certificate on the marquee when the certificate's limitation is the first thing a reviewer will test. The honest claim is: **a transparent, fully auditable BFS energy estimator** that achieves FoldX-neighbourhood Spearman ρ with zero force-field hand-tuning and zero learned parameters, and that makes every per-step energy attribution inspectable.

---

## 2. Correlation Performance

The numbers below are from the Phase 9-E MutCert run on the S2648 split (T4 lysozyme train n=97/val n=60, CI2 train n=27/val n=16, barnase train n=46/val n=31 after WT gating). Competitor figures are from published literature on various S2648 or S669 subsets and **cannot be directly compared** — they are shown for rough orientation only, not as evidence that MutCert matches or exceeds them. To make a valid comparison, all methods would need to run on the same split.

| Method | Family | Spearman ρ | Pearson r | n (val) | Notes |
|--------|--------|:----------:|:---------:|:-------:|-------|
| **MutCert v1** | T4 lysozyme | **−0.449** | −0.450 | 60 | This split |
| **MutCert v1** | Barnase | **−0.441** | −0.702 | 31 | This split, after WT gate |
| **MutCert v1** | CI2 | −0.194 | −0.368 | 16 | This split, n too small for reliable estimate |
| FoldX 5 | T4 lysozyme | ~−0.45† | ~−0.60† | ~100 | Different split |
| Rosetta ddg_mon | T4 lysozyme | ~−0.43† | ~−0.58† | ~100 | Different split |
| DDGun | T4 lysozyme | ~−0.40–0.44† | — | — | Parameter-free baseline; **must be run on same split** |
| ESM-1v | Multi-protein | ~−0.41† | ~−0.46† | ~500 | Different proteins, different split |
| MutPred2 | Multi-protein | ~−0.38† | ~−0.42† | ~500 | Different proteins, different split |

†Literature figures from Shu et al. 2020, Frenz et al. 2023, Meier et al. 2021 — not verified independently. Values may differ on this split.

> **Primary metric: Spearman ρ** — Pearson is inflated by outliers (barnase Pearson −0.702 vs Spearman −0.441 illustrates this; see `tribunal/opus-review-phase9d.md`). Always report both in any table.

> **DDGun is the most important missing baseline**: it is also parameter-free, sequence-and-structure-based, and antisymmetric. Its absence makes the "zero parameters" argument look weaker. It must be run on the same split before any submission.

---

## 3. Capability Comparison

| Capability | MutCert v1 | FoldX | Rosetta | DDGun | ESM-1v | MutPred2 |
|------------|:----------:|:-----:|:-------:|:-----:|:------:|:--------:|
| Point ΔΔG estimate | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Per-residue energy attribution | **✓** | partial | partial | ✗ | ✗ | ✗ |
| Glass-box BFS trace (ordered, per-step) | **✓** | ✗ | ✗ | ✗ | ✗ | ✗ |
| Running interval [lo, hi] | **✓** | ✗ | ✗ | ✗ | ✗ | ✗ |
| Interval narrows monotonically (by construction) | **✓** | ✗ | ✗ | ✗ | ✗ | ✗ |
| Convergence record stored and auditable | **✓** | ✗ | ✗ | ✗ | ✗ | ✗ |
| Zero force-field hand-tuning | **✓** | ✗ | ✗ | ✓ | ✓ | ✓ |
| Zero learned parameters | **✓** | ✗ | ✗ | ✓ | ✗ | ✗ |
| Sub-second CPU inference | **✓** | ✗ | ✗ | ✓ | partial | ✓ |
| Certified interval contains true ΔΔG | **not guaranteed** | — | — | — | — | — |

**Note on "monotone narrowing"**: this is guaranteed **by construction** because the uncertainty term ε₀/√k is monotone decreasing. It is not an empirical claim about accuracy. The interval reliably narrows onto the running MJ sum, which may itself be far from the true ΔΔG (see the L121A case below).

---

## 4. Known Limitation: Precision Without Accuracy

The L121A trace in `data/trace_MUT_P00720_A_121_ALA.html` is the canonical example:

- Certificate ISSUED (width = 0.49 kcal/mol ≤ 2.0, C1–C3 all satisfied)
- MutCert interval: approximately [−2.47, −1.98] kcal/mol
- Experimental ΔΔG: approximately +22.8 kcal/mol (strongly destabilising)
- Interval is ~25 kcal/mol away from experiment

The MJ potential assigns favourable contact-energy change for alanine at this buried hydrophobic position; the ratchet faithfully converges onto that wrong answer. This is not a bug in the convergence machinery — it is the expected behaviour when the underlying potential makes a large error.

The honest way to present this: **Spearman ρ ~0.45 means the ranking is moderately correct; individual predictions can be badly wrong**. The certificate tracks convergence of the estimation procedure, not accuracy of the result. A preprint that leads with the certificate without this caveat will not survive review.

---

## 5. Competitor Deep-Dive

### 5.1 DDGun (Montanucci et al., 2019; Pancotti et al., 2022)

**Approach**: Sequence- and structure-based, antisymmetric ΔΔG prediction, no machine learning, no hand-tuned force field.  
**Why it matters most**: DDGun is the most direct comparator for MutCert. Both are parameter-free, both use structural information, both avoid learning. The paper cannot claim "zero parameters" as a distinguishing feature without running DDGun on the same split.  
**Action item**: run DDGun on the same S2648 T4/barnase val split before submission.

### 5.2 FoldX 5 (Schymkowitz et al., 2005)

**Approach**: Physics-based energy function, ~30 hand-tuned terms, full side-chain repacking.  
**Strengths**: Best-in-class on buried mutations; widely cited; ~1 s per mutation.  
**Weaknesses**: Hand-tuned weights (not zero-parameter); black-box energy attribution; no convergence tracking.  
**MutCert difference**: comparable Spearman at zero tuning cost; full per-residue attribution; convergence record.

### 5.3 Rosetta ddg_monomer (Kellogg et al., 2011; Park et al., 2016)

**Approach**: Monte Carlo with REF2015 energy function.  
**Strengths**: High accuracy on buried mutations; large ecosystem.  
**Weaknesses**: Slow (minutes to hours per mutation for robust averaging); not auditable per-residue.  
**MutCert difference**: sub-second inference; fully deterministic; all attribution stored.

### 5.4 ESM-1v / ESM-2 (Meier et al., 2021; Lin et al., 2023)

**Approach**: Masked language model on ~250 M sequences; ΔΔG ≈ −Δlog P.  
**Weaknesses**: No structural mechanism; no attribution; point estimate only.  
**MutCert difference**: MutCert attributes every kcal/mol to a specific residue at a known hop distance. ESM cannot answer "which contact contributed most?"

### 5.5 MutPred2 (Li et al., 2019)

**Approach**: Gradient boosting on 1400+ features; optimised for pathogenicity, not ΔΔG.  
**Weaknesses**: Poor quantitative ΔΔG; feature importance ≠ mechanism; no convergence.  
**MutCert difference**: MutCert annotations come from deterministic BFS traversal, not feature correlation.

---

## 6. The Honest Smaller Claim

The story that survives review is narrower than "certified measurement protocol" — but it is real:

> MutCert is a **transparent, parameter-free BFS energy estimator** built on the Miyazawa–Jernigan contact potential with a reference-state correction. It achieves Spearman ρ ~0.44–0.45 on T4 lysozyme and barnase against the S2648 benchmark — in the same neighbourhood as FoldX and Rosetta — using zero force-field hand-tuning. Its distinguishing feature is that every step of the energy calculation is recorded and auditable: the glass-box trace shows which residues contributed, in what order, at what hop distance, and by how much.

> The convergence-tracking heuristic (the ratchet interval) is a useful diagnostic that makes the convergence behaviour of the BFS sum visible. It certifies precision, not accuracy. Improving the potential function (Phase 2) is the path to accuracy.

This claim:
- Is true
- Is differentiable from FoldX (attribution) and DDGun (trace)
- Does not require defending L121A as a "certified" result
- Sets up Phase 2 improvements naturally

---

*Generated: Phase 11 — Preprint Prep*  
*Keywords: differentiation, ddgun, foldx, rosetta, esm, mutpred2, glass-box-trace, transparency, spearman-449, spearman-441, l121a-limitation, precision-not-accuracy, honest-claim, phase11*
