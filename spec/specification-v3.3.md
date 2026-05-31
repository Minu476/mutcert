# MutCert — Convergence-Certified Mutation Stability Prediction

_version: 3.3 | status: Final | owner: Nasser Towfigh | updated: 2026-05-05_

**Project short name:** MutCert
**Full title:** Predicting thermodynamic stability change ($\Delta\Delta G$) caused by single amino acid mutations using graph-native Rich Learning agents with empirically validated convergence certificates.

A research project applying the **Rich Learning paradigm** to a specific, publishable, and patent-defensible problem in computational biology: **predicting the thermodynamic stability change ($\Delta\Delta G$) caused by a single amino acid mutation — with an empirically bounded convergence interval and a convergence certificate, not a point estimate.**

This is not a re-implementation of AlphaFold. AlphaFold predicts a folded structure from a sequence. This framework takes a *known* structure, applies a mutation, and predicts *by how much* stability changes — and crucially, *why*, with a full agent trace. AlphaFold DB is used as a cheap structural prior to bootstrap the graph, not as a competitor to beat.

---

## Research Paradigm

The Rich Learning paradigm applied here rests on four pillars:

**1. Natural graph representation.**
Proteins are inherently three-dimensional geometric graphs. Amino acids act as nodes; chemical bonds and spatial interactions (H-bonds, disulfide bridges, hydrophobic contacts, Van der Waals proximity, electrostatic charges) act as typed edges. Mutation effect prediction is a *local graph perturbation problem*: swap one node, let agents propagate the energy-change signal through the surrounding edge network, and read out the $\Delta\Delta G$ bound when the signal converges.

**2. Modular knowledge isolation — per-family skill modules.**
Each protein family (kinases, GPCRs, globins, etc.) is stored in its own dedicated, isolated Neo4j database within the DBMS. When a new family arrives, it is grafted by executing a `CREATE DATABASE {familyId}` command. Prior families are never retrained, queried, or modified during this process. This is a **structural, hardware-level isolation guarantee** — cross-database writes are architecturally impossible in Neo4j DBMS. This is distinct from continual learning "zero forgetting" claims: MutCert does not learn a shared model across families. Each family database is a permanently fixed, independently calibrated knowledge module.

**3. Empirically bounded convergence intervals — the core differentiator.**
Neural networks output a point estimate. MutCert agents output a **bounded convergence interval**: the system reports when the interval width has narrowed below a threshold and issues a convergence certificate. The interval bounds are heuristic convergence criteria, calibrated empirically per family. The scientific claim is empirical and falsifiable: $\ge 80\%$ coverage on held-out mutations, validated by a calibration (reliability) curve.

**4. Glass-box transparency.**
Every agent step — which residue interaction was re-evaluated, what the energy delta was, whether the edge weight change propagated or was absorbed — is recorded in the FSDE run graph. A researcher can replay the exact causal chain from mutation to final $\Delta\Delta G$ interval.

---

## Problem Statement

Single-point mutations in proteins are the primary mechanism behind genetic disease, drug resistance in pathogens, and therapeutic antibody engineering. Current methods fall into two categories:

- **Neural network models** (MutPred, ThermoNet, ESM variants): high average accuracy on benchmarks, but produce point estimates with no uncertainty bound, no convergence guarantee, and no causal explanation.
- **Physics-based simulators** (FoldX, Rosetta ddg_monomer): more interpretable but slow, and their energy functions are hand-tuned. FoldX achieves ~70–75% on S2648 at a ±1.0 kcal/mol window.

MutCert's claim is narrow and testable: for protein families where at least one experimental $\Delta\Delta G$ dataset exists, a Rich Learning graph agent can produce a bounded $\Delta\Delta G$ interval that contains the experimental value at $\ge 80\%$ frequency, with the interval width converging monotonically as more edges are evaluated — **certified and empirically validated, not assumed**.

Phase 1 threshold (±2.0 kcal/mol) is intentionally more lenient than FoldX's ±1.0 kcal/mol — the goal is to establish the certificate mechanism first, then tighten in phase 2.

---

## Objectives

- Build a data pipeline: UniProt sequence $\rightarrow$ AlphaFold DB structure $\rightarrow$ Neo4j DBMS (structural prior); overlay experimental $\Delta\Delta G$ values from ThermoMutDB/S2648 (primary) and BenchStab/S669 (secondary).
- Represent each protein family as an isolated Neo4j database; graft new families by creating new databases without touching existing ones.
- Maintain a central `run_registry` database to hold `ConvergenceCertificate` nodes and cross-family metadata.
- Implement the mutation agent: swap one residue node, re-evaluate all incident edges per the edge topology rules table, propagate energy-change signal via MJ contact potentials using BFS with distance decay, update the $\Delta\Delta G$ interval.
- Implement the convergence checker: track interval width per mutation step; issue convergence certificate when width is non-increasing for $\ge 5$ consecutive steps, width $\le 2.0$ kcal/mol, and all edges in an 8Å shell are evaluated.
- Validate bounded intervals: produce a **calibration (reliability) curve** — coverage % vs interval width threshold — and compare against FoldX and MJ-direct baselines.
- Record every agent run in the FSDE run graph; apply RunStep retention policy after certificate emission.

---

## Phase 1 Target Protein Families

| Family | UniProt | CIF source | ΔΔG records (ThermoMutDB) |
|---|---|---|---|
| **T4 lysozyme** | P00720 | RCSB 2LZM | 574 single-point |
| **CI2** | P01053 | AlphaFold EBI v6 | TBD |
| **barnase** | P00648 | AlphaFold EBI v6 | TBD |

---

## Datasets

| Dataset | Role | File on WD-Black | Records |
|---|---|---|---|
| **ThermoMutDB (single-point)** | Primary ΔΔG ground truth | `s2648/s2648.csv` | 8,806 |
| **BenchStab 2024** | Secondary / independent benchmark | `s669/s669.csv` | 289 |
| **ThermoMutDB (full)** | Full curated archive | `thermomutdb/thermomutdb.csv` | 13,337 |
| ProThermDB | Tertiary — CSV exports only if available | — | — |

Attribution: Potapov et al. 2009 (S2648 original); Velecký et al. 2024, *Bioinformatics* 40(9) (BenchStab); Nalisnick & Smyth (ThermoMutDB).

---

## Convergence Certificate Mechanism

**Interval initialisation (per mutation):**
The initial interval $[lo_0, hi_0]$ for a new mutation is set to:
$$[lo_0, hi_0] = [\Delta\Delta G_{MJ} - \varepsilon_0,\ \Delta\Delta G_{MJ} + \varepsilon_0]$$
where $\Delta\Delta G_{MJ}$ is the raw Miyazawa-Jernigan delta computed at the mutation site (sum of $e(mutant, neighbour) - e(wildtype, neighbour)$ over all edges incident to the mutated residue), and $\varepsilon_0$ is a per-family calibration parameter initialized from the training set variance:
$$\varepsilon_0 = \sigma(\Delta\Delta G_{train,\ family})$$

**Interval narrowing (per agent step):**

$$lo_{new} = \max(lo_{prev},\ \hat{\Delta\Delta G}_k - \varepsilon_{step}(k))$$
$$hi_{new} = \min(hi_{prev},\ \hat{\Delta\Delta G}_k + \varepsilon_{step}(k))$$

where $\hat{\Delta\Delta G}_k$ is the running sum of edge-weight deltas after $k$ edges evaluated, and:
$$\varepsilon_{step}(k) = \varepsilon_0 / \sqrt{k}$$

This guarantees monotone narrowing.

**Certificate emission:**
A `ConvergenceCertificate` node is written to the `run_registry` database when all three conditions hold:

1. $width(k) = hi(k) - lo(k) \le 2.0$ kcal/mol
2. $width(k) \le width(k-1) \le \ldots \le width(k-4)$ (non-increasing for 5 consecutive steps)
3. All edges in the K-nearest-neighbour shell (radius $\le 8$ Å) of the mutated residue have been evaluated at least once.

**Calibration note:**
$\varepsilon_0$ is learned per family from the 80% training split. The 2.0 kcal/mol threshold is a phase 1 starting point; the target is tightened to 1.0 kcal/mol in phase 2 by refining $\varepsilon_0$ and MJ weights per family.

---

## Validation Methodology

The primary validation output is a **calibration (reliability) curve**: for a sweep of interval width thresholds $w \in \{0.5, 1.0, 1.5, 2.0, 2.5, 3.0\}$ kcal/mol, compute the fraction of held-out mutations for which the issued interval contains the experimental $\Delta\Delta G$ value. A well-calibrated system produces a curve that tracks the diagonal.

**Required comparisons:**

| Method | Expected coverage at ±2.0 kcal/mol |
|---|---|
| Mean predictor baseline (family mean ΔΔG) | ~50% |
| MJ direct estimate (no BFS propagation) | ~60–65% |
| FoldX | ~70–75% |
| **MutCert (target)** | **≥ 80%** |

**Holdout split:** random 80/20 stratified by mutation severity ($|\Delta\Delta G| < 1$ / $1–3$ / $> 3$ kcal/mol). Split indices fixed in `data/s2648_split.json` before any validation run.

---

## Graph Edge Creation Rules

When loading a `.mmCIF` structure into Neo4j, the following rules determine which edge type to create between residue pair (i, j). Rules applied in priority order; a pair may have more than one edge type.

| Edge type | Neo4j label | Creation criterion |
|---|---|---|
| Peptide bond | `[:PEPTIDE]` | Consecutive residues in same chain — always created |
| H-bond | `[:H_BOND]` | N/O heavy-atom distance ≤ 3.5 Å AND donor-H-acceptor angle ≥ 120°; N, O, S are donor/acceptor atoms |
| Disulfide bridge | `[:DISULFIDE]` | Both residues Cys; Sγ–Sγ distance ≤ 2.1 Å; Cβ-Sγ-Sγ-Cβ dihedral $|\chi| \approx 90°\ (\pm 30°)$ |
| Hydrophobic contact | `[:HYDROPHOBIC_CONTACT]` | Both residues hydrophobic (Ala, Val, Leu, Ile, Met, Phe, Trp, Pro, Tyr); Cβ–Cβ distance ≤ 8.0 Å |
| Electrostatic pair | `[:ELECTROSTATIC]` | At least one residue charged (Asp−, Glu−, Lys+, Arg+, His±); side-chain centroid distance ≤ 10.0 Å |
| Van der Waals | `[:VAN_DER_WAALS]` | Any residue pair; Cα–Cα distance 3.0–8.0 Å; not already linked by a stronger interaction |

All distances from heavy atoms using AlphaFold/RCSB `.mmCIF` coordinates. H atoms absent in AlphaFold structures; H-bond criterion uses N/O heavy-atom distance only.

---

## Residue Swap and Edge Re-evaluation Rules

Swapping a residue changes edge topology, not just node properties. The mutation agent applies the following rules per residue-pair transition:

| Mutation type | Edge action |
|---|---|
| Any → Pro | Delete all `[:H_BOND]` edges where the mutated residue is the NH donor (Pro has no backbone NH) |
| Any → Cys or Cys → any | Re-evaluate all `[:DISULFIDE]` edges; create if Cβ–Cβ ≤ 7.5 Å and dihedral ≈ 90°; delete existing if source is no longer Cys |
| Any → charged (Asp, Glu, Lys, Arg, His) or reverse | Create/delete `[:ELECTROSTATIC]` edges to all residues within 10 Å |
| Any → Gly | Delete steric-clash constraints on alpha carbon; re-evaluate `[:VAN_DER_WAALS]` edges in 5 Å shell |
| Any → any (general) | Delete all `[:HYDROPHOBIC_CONTACT]` edges incident to the mutated node; recreate based on MJ hydrophobicity classification of the new residue |

---

## Energy Propagation Rule

**Contact potential: Miyazawa-Jernigan (MJ) 20×20 matrix.**
Computed from experimental PDB structures (not hand-tuned); residue-pair-specific; O(1) table lookup per edge. Non-MJ edge types influence *topology* only — their energy contribution to $\hat{\Delta\Delta G}$ is approximated via MJ lookup for dimensional consistency. Phase 2 may introduce additive Coulomb/SASA corrections with explicit unit conversion.

**Propagation algorithm (BFS with distance-weighted decay):**

1. Start at the mutated residue node M.
2. Compute $\Delta\Delta G_{local} = e(mutant, neighbour) - e(wildtype, neighbour)$ for each edge incident to M.
3. Push affected neighbours into BFS queue Q.
4. For each node N at graph distance $d$ from M:
   - Recompute all edges incident to N.
   - Apply decay: $\Delta e_{eff} = \Delta e_{raw} \times \exp(-d / \lambda)$, $\lambda = 2.5$ hops $\approx 8$ Å.
   - Accumulate into $\hat{\Delta\Delta G}$.
   - Push N's unvisited neighbours if $\Delta e_{eff} > \varepsilon_{cutoff} = 0.01$ kcal/mol.
5. Stop when Q is empty or all nodes beyond the 8 Å shell are below $\varepsilon_{cutoff}$.

**Scope and limitations:**
BFS propagation is a local approximation. Long-range allosteric effects beyond the 8 Å shell are not captured in phase 1. This is an explicit design constraint, not an oversight — the phase 1 goal is to establish the certificate mechanism on the local shell, then extend in phase 2.

---

## RunStep Retention Policy

After a `ConvergenceCertificate` is issued for a `MutationRun`, the detailed `RunStep` nodes are compressed to reduce Neo4j storage:

- **Keep permanently:** first `RunStep`, last `RunStep`, any `RunStep` flagged with a rejection code (e.g. `STERIC_CLASH`, `INTERVAL_INVERSION`), and the `ConvergenceCertificate` node itself.
- **Archive:** all intermediate `RunStep` nodes — strip detailed edge properties, retain only `{stepIndex, energyDelta_kcal, graphDistance, accepted}`.
- **Trigger:** compression runs automatically after certificate emission, within the same write transaction.
- **Replay guarantee:** the first + last + flagged steps are sufficient to reconstruct the causal summary; full replay of intermediate steps is not guaranteed post-compression.

---

## Architecture

### Neo4j DBMS Layout

```text
neo4j DBMS (Desktop Enterprise)
├── db: t4_lysozyme      ← Family 1 (Isolated nodes/edges)
├── db: ci2              ← Family 2 (Isolated nodes/edges)
├── db: barnase          ← Family 3 (Isolated nodes/edges)
├── db: run_registry     ← Cross-family index, ConvergenceCertificates, FSDE states
└── db: system           ← Built-in
```

### Memory Tuning (`neo4j.conf`)

Required before Phase 1 execution (48GB RAM workstation):

```properties
server.memory.heap.initial_size=4g
server.memory.heap.max_size=8g
server.memory.pagecache.size=12g
```

### Neo4j Schema (outline)

```cypher
(:Residue {id, chain, seqPos, aminoAcid, phi, psi, charge, hydrophobicity})
  -[:PEPTIDE {length_A}]->
  -[:H_BOND {energy_kcal, donor, acceptor}]->
  -[:HYDROPHOBIC_CONTACT {mj_energy_kcal}]->
  -[:ELECTROSTATIC {charge_product, distance_A}]->
  -[:DISULFIDE {cbeta_dist_A, dihedral_deg}]->
  -[:VAN_DER_WAALS {mj_energy_kcal}]->

(:ProteinFamily {uniprotId, pfamId, name, importedAt, epsilon0, mjWeights})
  -[:HAS_RESIDUE]->(:Residue)

(:MutationRun {id, wildType, mutant, seqPos, familyId, timestamp})
  -[:PRODUCED]->(:DeltaDeltaGInterval {lo, hi, width, certified, experimentalValue})
  -[:RUN_STEP {stepIndex}]->(:RunStep {action, edgeType, energyDelta_kcal, accepted, graphDistance})

(:ConvergenceCertificate {runId, width_kcal, stepsAtConvergence, shellRadiusA, timestamp})
```

### Transaction Isolation

Each `MutationRun` executes in its own Neo4j write transaction scoped strictly to its `RunStep` subgraph. Reads from the parent family database are performed in a separate read transaction before the write transaction opens. Concurrent agents on the same family do not block each other.

---

## Tech Stack

### Allowed

- **.NET 10 (C#)** — primary runtime; parser, mutation agent, Neo4j integration.
- **System.Threading.Channels / TPL Dataflow** — concurrent mutation agent runner.
- **Neo4j 5.x Enterprise** — persistent graph store, per-family databases.
- **FSDE Engine** — agent orchestration, session memory, run recording.
- **BioCif / custom parser** — streaming `.mmCIF` tokenisation.
- **Miyazawa-Jernigan 20×20 matrix** — embedded static lookup table in C# (Miyazawa & Jernigan 1996, Table 3).
- **Python scripts** (`scratch/`) — ThermoMutDB/BenchStab CSV preprocessing and EDA.
- **AlphaFold EBI REST API** — structural prior fetch.
- **UniProt REST API** — FASTA sequence and family annotation.
- **RCSB PDB REST API** — experimental structure fallback.

### Forbidden

- Black-box ML model calls during agent reasoning.
- Storing raw `.cif` / `.pdb` bytes in Neo4j.
- Modifying existing family database nodes when grafting a new family.
- Point-estimate ΔΔG output without an accompanying interval and certificate.

---

## Constraints

- All internal energy units: kcal/mol.
- `.mmCIF` parsing must be streaming (files can exceed 100 MB).
- Grafting a new family must not issue any Cypher `MATCH ... SET` or `MERGE ... SET` on nodes in prior family databases. Verified by node-count assertion before and after every graft.
- Convergence certificate emitted only when all three conditions hold (width, monotone, shell complete).
- API rate limits: RCSB 10 req/s, UniProt 1 req/s, AlphaFold EBI 10 req/s.
- Agent loop must be pausable and resumable via FSDE session commands.
- Validation split indices fixed in `data/s2648_split.json` before any validation run.

---

## Timeline

| Date | Milestone | Status |
|---|---|---|
| 2026-05-04 | Specification v3.0–v3.2 finalised; FSDE running | ✅ Done |
| 2026-05-05 | Data pipeline complete — FASTAs, CIFs, ThermoMutDB, BenchStab on WD-Black | ✅ Done |
| 2026-05-05 | Spec v3.3 finalised (Phase 0 complete) | ✅ Done |
| 2026-05-11 | `.NET 10` `Protein.Engine` scaffolded; Neo4j databases created; MJ matrix embedded | — |
| 2026-05-18 | `StructureLoader` complete; T4 lysozyme loaded into Neo4j; edge detection working | — |
| 2026-05-25 | BrainCore guard conditions + residue-swap edge re-evaluation rules implemented | — |
| 2026-06-01 | `MutationAgentPool` + `EnergySignalPropagator` running on T4 lysozyme | — |
| 2026-06-08 | `ConvergenceChecker` + certificate emission working; calibration curve plotted | — |
| 2026-06-15 | Coverage ≥ 80% on T4 lysozyme held-out set; reliability diagram validated | — |
| 2026-06-22 | CI2 + barnase grafted and validated; zero-forgetting verified | — |
| 2026-06-29 | Draft paper outline + patent claim structure | — |

---

## Differentiation Summary

| Capability | Neural nets (MutPred, ESM) | FoldX / Rosetta | **MutCert** |
|---|---|---|---|
| Bounded ΔΔG interval with certificate | No | No | **Yes** |
| Modular knowledge isolation (per family) | No | N/A | **Yes** |
| Calibration curve + reliability diagram | No | No | **Yes** |
| Causal step trace, fully replayable | No | Partial | **Yes** |
| Competes with AlphaFold | Yes | No | **No — uses it as prior** |
| C# / Neo4j native | No | No | **Yes** |
| CPU-native inference | No | Partial | **Yes** |
| Phase 1 coverage baseline | ~82% at ±2.0 (ESM) | ~73% at ±1.0 | **≥ 80% at ±2.0 (target)** |

---

## Security & Compliance

- No API keys needed for phase 1 (all APIs are public/unauthenticated).
- Downloaded files validated (size + format check) before parsing.
- No patient data or proprietary sequences — CC0/CC-BY licensed data only.
- Required attributions: Jumper et al. 2021 (AlphaFold); Miyazawa & Jernigan 1996 (MJ matrix); Velecký et al. 2024 (BenchStab); ThermoMutDB curators.
