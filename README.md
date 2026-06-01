# MutCert

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20499591.svg)](https://doi.org/10.5281/zenodo.20499591)

**A transparent, parameter-free graph estimator for protein stability changes (ΔΔG) of single-point mutations.**

MutCert loads a protein structure as a residue **contact graph**, and when a residue is mutated it walks the graph outward from the mutation site, accumulating Miyazawa–Jernigan pairwise contact energies hop-by-hop with distance decay. Every step is recorded, so the prediction is fully auditable: you can see exactly which neighbouring residue contributed which fraction of a kcal/mol, and in what order.

It has **zero learned parameters and zero force-field hand-tuning** — only the published Miyazawa–Jernigan 20×20 contact-energy table (1996) plus a Swiss-Prot reference-state correction, and three closed-form calibration constants (ε₀ and an OLS slope/intercept) fit per family.

---

## ⚠️ What the convergence certificate does and does not mean

MutCert emits a `ConvergenceCertificate` when its running interval `[lo, hi]` has narrowed below a width threshold and the contact shell is exhausted. **This certifies that the Miyazawa–Jernigan energy sum has converged — it does *not* certify that the prediction is correct.**

The interval is a **convergence interval of the MJ sum, not a calibrated uncertainty interval on the true ΔΔG.** A clear illustration ships with the repo: the **L121A** mutation in T4 lysozyme issues a certificate with interval width **0.49 kcal/mol**, while the experimental ΔΔG lies roughly **25 kcal/mol outside** that interval (see `data/trace_MUT_P00720_A_121_ALA.html`). The ratchet faithfully converges onto a value that the MJ potential simply gets wrong. Treat the certificate as a statement about numerical convergence and auditability, never as a confidence bound on accuracy.

The reported **convergence rate (98.9%, 932 of 942 runs)** measures only how often the MJ sum converged, not how often the prediction was accurate. This count covers all mutation runs across all three families (train+val splits plus exploratory runs), not just the 107 validation mutations.

---

## How it works

1. **Structure → contact graph.** AlphaFold (or RCSB) mmCIF structures are parsed and loaded into Neo4j as residue nodes with typed contact edges (peptide, hydrogen-bond, hydrophobic, electrostatic, van der Waals, disulfide). One database per protein family keeps families fully isolated.
2. **BFS energy propagation.** From the mutation site, a breadth-first walk visits contacting residues in order of hop distance, computing each step's MJ contact-energy change (with a per-residue Swiss-Prot reference-state correction) and an exponential hop-distance decay.
3. **Running interval [lo, hi].** A monotone "ratchet" narrows the interval as more contacts are evaluated; a certificate is issued once the width is below threshold and the contact shell is exhausted.
4. **Glass-box trace.** Every BFS step — residue visited, hop distance, energy contribution, and the running interval — is logged and rendered as a self-contained, inspectable HTML report.

A per-family OLS recalibration (`ΔΔG ≈ α + β·MJ_direct`) is fit on the training split to correct the systematic scale/sign offset of raw MJ sums; see *Limitations* for what this does and does not fix.

---

## Results

Validated on the **S2648** benchmark (Potapov et al., 2009), stratified 80/20 split, three protein families. The split is frozen in `data/s2648_split.json` so results reproduce exactly.

| Family | Val n | MutCert Spearman \|ρ\| | MutCert Pearson \|r\| | DDGun-seq Spearman \|ρ\| (same split) |
|---|---:|---:|---:|---:|
| T4 lysozyme | 60 | 0.449 | 0.450 | 0.426 |
| Barnase | 31 | 0.441 | 0.702† | 0.731 |
| CI2 | 16 | 0.194‡ | 0.368‡ | 0.635 |

†Barnase Pearson is inflated by a few high-magnitude outliers; rank correlation (Spearman) is the honest measure.  
‡CI2 (n=16) is underpowered — 95% CI for ρ spans roughly [−0.3, +0.6], indistinguishable from zero.

**Spearman ρ is the primary metric.** Pearson is reported alongside it for transparency, but is not the headline because barnase's Pearson is inflated by a few high-magnitude outliers; the rank correlation is the honest measure of predictive skill.

**Honest reading:**
- On **T4 lysozyme**, MutCert reaches DDGun-seq's neighbourhood (0.449 vs 0.426 at n=60 — statistically indistinguishable).
- On **barnase** (the best-powered comparison), **DDGun-seq is clearly stronger** (0.731 vs 0.441).
- **CI2 is underpowered.** At n=16, a Spearman of 0.194 is not statistically distinguishable from zero (its 95% CI spans roughly −0.3 to +0.6), so it should be read as *inconclusive*, not as a result in either direction.

MutCert is **not** claimed to beat other methods on accuracy. Its contribution is **transparency**: a deterministic, parameter-free, per-step-auditable estimate, with a convergence interval and a complete causal trace that point methods (FoldX, DDGun, ESM-1v) do not provide.

---

## Requirements

- **.NET 10 SDK**
- **Neo4j 5.x Enterprise Edition** — the per-family model uses **multiple named databases**, which requires Enterprise (Community edition supports only a single user database)
- Docker (optional, for containerized setup)

Neo4j connection is read from environment variables:
- `NEO4J_URI` (e.g., `bolt://localhost:7687`)
- `NEO4J_USERNAME` (default: `neo4j`)
- `NEO4J_PASSWORD` (default: `mutcert`)

**No credentials are committed to this repository.** See `.env.example` for setup.

## Install & build

### Docker (recommended)

```bash
git clone https://github.com/nassertowfigh/mutcert
cd mutcert
cp .env.example .env   # defaults work out of the box
docker compose up -d neo4j
# Wait ~30s for Neo4j to be ready
docker compose run --rm mutcert-init   # one-time: graft protein families
```

### Native (.NET)

```bash
git clone https://github.com/nassertowfigh/mutcert
cd mutcert
dotnet build MutCert.sln

# Set environment variables
export NEO4J_URI=bolt://localhost:7687
export NEO4J_USERNAME=neo4j
export NEO4J_PASSWORD=<your-password>

# Graft protein families (one-time setup)
dotnet run --project src/Protein.Engine -- graft P00720 t4-lysozyme PF00959 data/cif/t4_lysozyme_P00720_2LZM.cif
dotnet run --project src/Protein.Engine -- graft P01053 ci2 PF00014 data/cif/ci2_P01053.cif
dotnet run --project src/Protein.Engine -- graft P00648 barnase PF00211 data/cif/barnase_P00648.cif
dotnet run --project src/Protein.Engine -- import-mutations
```

## Usage

```bash
# Validate all three families end-to-end (regenerates data/output/validation_report.html)
dotnet run --project src/Protein.Engine -- validate-all

# Validate one family
dotnet run --project src/Protein.Engine -- validate-family <family>

# Score a single mutation
dotnet run --project src/Protein.Engine -- run-mutation <MUT_ID>

# Render the glass-box trace for one mutation
dotnet run --project src/Protein.Engine -- trace-mutation <MUT_ID> [outputPath]

# Verify deterministic per-family database isolation
dotnet run --project src/Protein.Engine -- verify-isolation
```

Full subcommand list: `graft`, `import-mutations`, `run-mutation`, `batch-run`, `validate-family`, `validate-all`, `trace-mutation`, `replay`, `verify`, `verify-isolation`, `whereami`.

## Outputs

- `data/output/validation_report.html` — per-family metrics, calibration curves, and an embedded 3D structure viewer.
- `data/output/trace_<MUT_ID>.html` — the full glass-box causal trace for a single mutation, including a convergence chart and per-step energy breakdown.

---

## Data & provenance

| Asset | Source | License / status | Attribution |
|---|---|---|---|
| Protein structures (CIF) | AlphaFold DB / RCSB (2LZM) | AlphaFold: **CC-BY 4.0**; RCSB coordinates: public | DeepMind & EMBL-EBI (AlphaFold); RCSB PDB (2LZM) |
| S2648 mutations | Potapov et al., 2009 | **Must be obtained by user** (see `data/s2648/MANUAL_DOWNLOAD.txt`) | Potapov, Cohen & Schreiber, 2009 |
| S669 mutations | Pancotti et al., 2022 | **Must be obtained by user** (see `data/s669/MANUAL_DOWNLOAD.txt`) | Pancotti et al., 2022 |
| Miyazawa–Jernigan contact potential | Miyazawa & Jernigan, 1996, *J. Mol. Biol.* 256:623–644 | published table (data/fact) | Miyazawa & Jernigan, 1996 |
| DDGun-seq baseline | Montanucci, Fariselli et al. | DDGun 0.0.2 (patched for Python 3.13 / BioPython 1.87) | DDGun authors |

> S2648 and S669 CSV files are **not redistributed** in this repository due to unclear licensing terms. Users must download them from original sources (ThermoMutDB or paper supplementary materials). See `MANUAL_DOWNLOAD.txt` files for instructions and expected checksums.

## Reproducibility

- **Frozen train/val split**: `data/s2648_split.json` (stratified by severity tier, committed to repo).
- **Deterministic**: No stochastic elements; BFS order, OLS fit, and certificate logic are all deterministic, so a given run is bit-for-bit reproducible.
- **DDGun-seq comparison**: Produced on the *same* split via `scripts/ddgun_baseline.py` using DDGun 0.0.2 with ColabFold/MMseqs2 MSA profiles (1180–4589 sequences per family). Results saved in `data/ddgun_baseline_results.json` (on external storage, not committed — see data/ symlink to WD-Black).

## Limitations

- Predictive accuracy is **modest** (Spearman ~0.45 on the one experimentally-structured family) and below sequence/structure-learning methods such as ESM-based predictors.
- The OLS recalibration corrects scale and sign offset but **cannot improve rank correlation** (Spearman is invariant under monotonic transforms); higher accuracy would require a richer energy model (explicit H-bond and electrostatic terms), not more calibration.
- The convergence certificate is **decoupled from accuracy** (see the L121A caveat above).
- S2648 is a forward-mutation benchmark; the antisymmetry-testing, leakage-controlled **S669** would be a more rigorous evaluation in future work.
- **"Per-family isolation"** means each family lives in its own Neo4j database, so adding a family cannot alter another's results — a **deterministic isolation property**, not a continual-learning ("no catastrophic forgetting") claim. We avoid the term "zero forgetting" to prevent confusion with machine-learning terminology.

## Citing

If you use MutCert, please cite via the metadata in [`CITATION.cff`](CITATION.cff) and the archived Zenodo record (DOI: *to be added upon publication*).

## License

MIT License — see [LICENSE](LICENSE) for full text.
