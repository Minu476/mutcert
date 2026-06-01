# MutCert

A transparent, parameter-free protein stability change (ΔΔG) estimator built on the Miyazawa–Jernigan contact potential. Every step of the energy calculation is recorded and fully auditable via the glass-box BFS trace.

**Preprint status**: Phase 11 docs in `tribunal/`. Canonical validation results below.

---

## Canonical Validation Results (Phase 9-E, S2648 benchmark)

| Family | Train n | Val n | Spearman \|ρ\| | Pearson \|r\| |
|--------|:-------:|:-----:|:------------:|:-----------:|
| T4 lysozyme | 97 | 60 | **0.449** | 0.450 |
| Barnase | 46 | 31 | **0.441** | 0.702 |
| CI2 | 27 | 16 | 0.194† | 0.368† |

†CI2 (n=16) is underpowered — 95% CI for ρ spans roughly [−0.3, +0.6], indistinguishable from zero.

> **Primary metric: Spearman ρ** — less sensitive to outliers than Pearson (e.g., barnase |r|=0.702 vs |ρ|=0.441). Sign convention: MutCert uses negative correlation (larger MJ sum → more stable), opposite to S2648 experimental convention.

**Run-registry:** 942 total mutation runs across all 3 families (T4/CI2/barnase train+val, plus exploratory runs), 932 certificates issued (98.9% convergence rate). This covers the full dataset, not just the 107 val mutations.

---

## DDGun Baseline Comparison (Same S2648 Val Split)

| Family | n | DDGun-seq \|ρ\| | MutCert \|ρ\| | Verdict |
|--------|:-:|:-------------:|:------------:|:-------:|
| T4 lysozyme | 60 | 0.426 | **0.449** | Tied |
| CI2 | 16 | **0.635** | 0.194 | DDGun better |
| Barnase | 31 | **0.731** | 0.441 | DDGun better |

DDGun-seq outperforms MutCert on 2 of 3 families. DDGun benefits from evolutionary MSA signal (1180–4589 sequences per family via ColabFold/MMseqs2); MutCert uses only local MJ contact potentials with no sequence database dependency.

**MutCert's differentiator is auditability, not raw accuracy:** every per-step energy attribution is recorded and inspectable via glass-box trace. No other method (FoldX, DDGun, ESM, Rosetta) provides ordered, hop-by-hop energy decomposition.

Script: `scripts/ddgun_baseline.py` reproduces these numbers (requires `ddgun` 0.0.2, Python 3.13+, BioPython 1.87+).

---

## Known Limitation: Precision ≠ Accuracy

The L121A trace (`MUT_P00720_A_121_ALA`) is the canonical example:
- Certificate **ISSUED** (width 0.49 kcal/mol ≤ 2.0)
- MutCert interval: [−2.47, −1.98] kcal/mol
- Experimental ΔΔG: **+22.8 kcal/mol** (strongly destabilizing)
- **Interval is ~25 kcal/mol away from truth**

The MJ potential assigns favorable contact-energy change for this buried hydrophobic→alanine mutation; the convergence machinery faithfully narrows onto that wrong answer. **The certificate guarantees convergence of the MJ sum (precision), not accuracy of the prediction.** See `tribunal/phase11-cert-definition.md §6` for full discussion.

---

## Quick Start (Docker)

### Prerequisites
- Docker Desktop ≥ 4.x
- ~4 GB free RAM (for Neo4j)

### 1. Clone and configure

```bash
git clone https://github.com/nassertowfigh/mutcert
cd mutcert
cp .env.example .env   # defaults work out of the box
```

### 2. Start Neo4j

```bash
docker compose up -d neo4j
# Wait ~30 s for Neo4j to be ready (check: http://localhost:7474)
```

### 3. Build protein graphs (one-time setup, ~2 min)

```bash
docker compose run --rm mutcert-init
```

This grafts the three protein families into Neo4j and imports the S2648 mutations.

> **S2648 CSV**: The benchmark CSV (`data/s2648/s2648.csv`) must be obtained from
> Potapov et al. 2009 supplementary materials or ThermoMutDB. See
> `data/s2648/MANUAL_DOWNLOAD.txt`. Place it at `data/s2648/s2648.csv` before
> running init.

### 4. Run validation (reproduces canonical numbers, ~5–10 min)

```bash
docker compose run --rm mutcert validate-all
# Output: data/output/validation_report.html
```

### 5. Generate a glass-box trace

```bash
docker compose run --rm mutcert trace-mutation MUT_P00720_A_121_ALA
# Output: data/output/trace_MUT_P00720_A_121_ALA.html
```

---

## Running Without Docker

Requirements: .NET 10 SDK, Neo4j 5 Enterprise (multi-database).

```bash
# Set env vars
export NEO4J_URI=bolt://localhost:7687
export NEO4J_USERNAME=neo4j
export NEO4J_PASSWORD=<your-password>

# Build
dotnet build src/Protein.Engine/

# Graft protein families
dotnet run --project src/Protein.Engine -- graft P00720 t4_lysozyme PF00959 data/cif/t4_lysozyme_P00720_2LZM.cif
dotnet run --project src/Protein.Engine -- graft P01053 ci2 PF00014 data/cif/ci2_P01053.cif
dotnet run --project src/Protein.Engine -- graft P00648 barnase PF00211 data/cif/barnase_P00648.cif

# Import mutations
dotnet run --project src/Protein.Engine -- import-mutations

# Run all validation
dotnet run --project src/Protein.Engine -- validate-all
```

---

## CLI Reference

```
dotnet run -- graft <uniprotId> <familyName> <pfamId> <cifPath>
dotnet run -- verify
dotnet run -- verify-isolation
dotnet run -- import-mutations
dotnet run -- run-mutation <uniprotId> <seqPos> <mutant3Letter>
dotnet run -- batch-run <uniprotId> <familyDbName> [workerCount]
dotnet run -- whereami
dotnet run -- replay <runId>
dotnet run -- validate-family <familyName>
dotnet run -- validate-all
dotnet run -- zero-forgetting
dotnet run -- trace-mutation <mutationId> [outputPath]
```

---

## Key Algorithm

MutCert uses a BFS traversal of the protein contact graph, accumulating Miyazawa–Jernigan contact energy changes. The running interval narrows via a max/min ratchet:

$$\text{uncertainty}(k) = \frac{\varepsilon_0}{\sqrt{k}}$$

$$\text{lo}(k) = \max(\text{lo}(k-1),\ \Sigma_k - \text{uncertainty}(k))$$
$$\text{hi}(k) = \min(\text{hi}(k-1),\ \Sigma_k + \text{uncertainty}(k))$$

A certificate is issued when:
- **C1**: interval width ≤ 2.0 kcal/mol
- **C2**: convergence flag set by `ConvergenceSupervisor`
- **C3**: ≥ 5 BFS steps completed

**Critical:** The interval [lo, hi] bounds the *convergence of the MJ sum*, NOT the *error of the prediction*. It is not a confidence interval or calibrated uncertainty quantification. The certificate records precision (the MJ sum has converged), not accuracy (the MJ sum may be far from experimental ΔΔG). See L121A example above and `tribunal/phase11-cert-definition.md` for full discussion.

---

## Structure

```
src/Protein.Engine/        ← .NET 10 console app
  BrainCore.cs             ← orchestration
  EnergySignalPropagator.cs← max/min ratchet interval
  ConvergenceSupervisor.cs ← certificate issuance + Neo4j persistence
  GlassBoxTrace.cs         ← Phase 11: HTML trace generator
  ValidationPipeline.cs    ← S2648 benchmark runner
  StructureLoader.cs       ← mmCIF parser + graph construction
  MiyazawaJernigan.cs      ← MJ contact potential (20×20 matrix)
data/
  cif/                     ← mmCIF structure files (PDB)
  s2648/                   ← S2648 benchmark split
  s2648_split.json         ← train/val assignment
tribunal/
  phase11-methods.md       ← full methods specification
  phase11-cert-definition.md ← certificate formal definition
  phase11-differentiation.md ← comparison vs FoldX/Rosetta/DDGun/ESM
```

---

## Limitations

- **L121A (T4 Lysozyme)**: Certificate ISSUED with width 0.49 kcal/mol; experimental ΔΔG ~+22.8 kcal/mol. The MJ potential assigns favourable contact energy at this buried hydrophobic site; the ratchet converges faithfully onto the wrong answer.
- CI2 Spearman ρ −0.194 (n=16 val) is below significance threshold. Small sample.
- DDGun baseline on the same split is the most important missing comparison. Numbers in `tribunal/phase11-differentiation.md` are from different splits and cannot be directly compared.

---

## License

MIT — see `LICENSE`.

The S2648 benchmark dataset is from Potapov et al. (2009) and must be obtained separately. CIF files are from the RCSB Protein Data Bank (public domain).

---

## Citation

If you use MutCert in research, please cite:

```bibtex
@software{mutcert2026,
  author = {Towfigh, Nasser},
  title  = {MutCert: Transparent BFS Energy Estimator for Protein Stability},
  year   = {2026},
  url    = {https://github.com/nassertowfigh/mutcert}
}
```
