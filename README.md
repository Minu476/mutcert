# MutCert

A transparent, parameter-free protein stability change (ΔΔG) estimator built on the Miyazawa–Jernigan contact potential. Every step of the energy calculation is recorded and fully auditable via the glass-box BFS trace.

**Preprint status**: Phase 11 docs in `tribunal/`. Canonical validation results below.

---

## Canonical Validation Results (Phase 9-E, S2648 benchmark)

| Family | Train n | Val n | Spearman ρ | Pearson r |
|--------|:-------:|:-----:|:----------:|:---------:|
| T4 lysozyme | 97 | 60 | **−0.449** | −0.450 |
| Barnase | 46 | 31 | **−0.441** | −0.702 |
| CI2 | 27 | 16 | −0.194 | −0.368 |

Run-registry: 942 total runs, 932 certified (98.9% convergence rate).

> **Note on L121A**: The glass-box trace for `MUT_P00720_A_121_ALA` shows cert ISSUED (width 0.49 kcal/mol), but the experimental ΔΔG is ~+22.8 kcal/mol — ~25 kcal/mol outside the interval. The certificate tracks convergence of the MJ sum, not accuracy. See `tribunal/phase11-cert-definition.md §6`.

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

The certificate records **precision** (convergence of the MJ sum), not accuracy. See `tribunal/phase11-cert-definition.md` for full specification including the L121A limitation.

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
