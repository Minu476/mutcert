# MutCert — Methods Section: MJ Potential, ε₀ Calibration & BFS Propagation

**Phase 11 Preprint Prep | MutCert v1.0**  
**Keywords:** methods-section, miyazawa-jernigan, epsilon0, bfs-propagation, reference-state, ratchet-narrowing, contact-graph, phase11

---

## 1. Contact Graph Construction

### 1.1 Residue Node Definition

Each protein is represented as a graph $G = (V, E)$ where:
- Each node $v_i \in V$ corresponds to a single amino acid residue, identified by chain and sequence position
- Node attributes: residue type (3-letter code), chain ID, sequence position, coordinates (Cα)
- Graph is constructed from the AlphaFold mmCIF structure, parsed via `MmcifParser.cs`

### 1.2 Edge Types and Criteria

Edges encode physical interactions identified from the structure:

| Edge type | Label | Criterion |
|-----------|-------|-----------|
| Backbone bond | `PEPTIDE` | Sequential residues (|Δseq| = 1) |
| Hydrogen bond | `H_BOND` | Minimum N/O/S heavy-atom distance ≤ 3.5 Å between residues |
| Hydrophobic contact | `HYDROPHOBIC_CONTACT` | Cβ–Cβ distance ≤ 8.0 Å, both residues hydrophobic (MJ classification) |
| Electrostatic pair | `ELECTROSTATIC` | Sidechain centroid distance ≤ 10.0 Å, at least one residue charged |
| Van der Waals | `VAN_DER_WAALS` | Cα–Cα distance 3.0–8.0 Å, no other edge already registered between the pair |
| Disulfide | `DISULFIDE` | Sγ–Sγ distance ≤ 2.5 Å, Cβ–Sγ–Sγ–Cβ dihedral 60°–120° (CYS only) |

All edges are stored as directed pairs (both directions) in Neo4j. The graph is **frozen** at import time and never modified during mutation analysis — only the BFS traversal reads from it.

### 1.3 Sequence Offset Correction

AlphaFold CIF files use UniProt canonical numbering, which includes signal peptides and propeptides. The S2648 benchmark uses mature-protein (PDB/experimental) numbering. The offset is applied at CSV parse time:

| Protein | UniProt | Offset | Reason |
|---------|---------|:------:|--------|
| T4 lysozyme | P00720 | 0 | RCSB structure, no prefix |
| CI2 | P01053 | +1 | 1-residue propeptide in AlphaFold CIF |
| Barnase | P00648 | +47 | 47-residue signal+propeptide in AlphaFold CIF |

Empirically verified by matching 3 residues per protein between S2648 WT codes and CIF `label_seq_id` + `label_comp_id` (Phase 9-D).

---

## 2. Miyazawa–Jernigan (MJ) Contact Potential

### 2.1 Potential Choice Rationale

The Miyazawa–Jernigan (MJ) statistical contact potential (Miyazawa & Jernigan, 1996) is used for three reasons:

1. **Interpretability**: A 20×20 symmetric matrix — every entry is a single number attributable to a pair of amino acid types. There are no hidden layers, no learned features.
2. **Speed**: O(1) lookup per residue pair — enables sub-second BFS on full-protein graphs.
3. **Benchmark status**: MJ is the canonical baseline contact potential in the biophysics literature, making MutCert results directly comparable to classical methods.

The MJ matrix encodes the effective contact energy $e_{ij}$ (in RT units at 300 K) for amino acid types $i$ and $j$, derived from statistical analysis of observed contact frequencies in known protein structures relative to a random reference state.

### 2.2 Reference-State Correction

Naive MJ energies over-penalise mutations at well-connected buried sites because all contacts of the wild-type are counted in the sum. The reference-state correction subtracts the mean-field energy expected for any amino acid in an average protein environment:

$$\Delta G_{\text{ref}}(\alpha) = \sum_{j=1}^{20} f_j \cdot e_{\alpha j}$$

where $f_j$ is the background frequency of amino acid type $j$ in the Swiss-Prot database (non-redundant, all kingdoms).

The corrected MJ energy contribution of step $k$ (residue $v_k$ at hop $h_k$, with $n_h$ evaluated hop-1 neighbours) is:

$$\Delta E_k = e_{\text{mut}, v_k} - e_{\text{wt}, v_k} - n_h \cdot (\Delta G_{\text{ref}}(\text{mut}) - \Delta G_{\text{ref}}(\text{wt}))$$

The $n_h$ factor ensures the reference correction scales with the number of contacts already evaluated — avoiding double-counting.

**Implementation**: `MiyazawaJernigan.cs` — `GetReferenceEnergy(string residueType)`, pre-computed in static constructor from Swiss-Prot background frequencies.

### 2.3 MJ Matrix Source

MJ values from Table 3 of Miyazawa & Jernigan (1996), *Journal of Molecular Biology* 256:623–644. The matrix is embedded directly in `MiyazawaJernigan.cs` as a 20×20 `double[,]` array indexed by amino acid type. Units: kcal/mol (converted from RT at 300 K via RT = 0.5961 kcal/mol at 300 K).

---

## 3. BFS Propagation Rule

### 3.1 BFS Order

The agent performs breadth-first search starting from the mutation site $v_m$, visiting residues in order of increasing hop distance:

1. **Hop 0**: The mutation site itself ($v_m$) — MJ self-energy change applied
2. **Hop 1**: All residues with an edge to $v_m$ in the contact graph
3. **Hop k**: All residues at graph distance $k$ from $v_m$, ordered by edge type priority within each shell

Within each hop distance, residues are visited in **deterministic graph order** (Neo4j node creation order), ensuring reproducibility.

### 3.2 Interval Narrowing Rule

The interval is updated by a **max/min ratchet** (`EnergySignalPropagator.cs`):

$$\text{uncertainty}(k) = \varepsilon_0 / \sqrt{k}$$
$$\text{lo}(k) = \max\bigl(\text{lo}(k-1),\; \Sigma_k - \text{uncertainty}(k)\bigr)$$
$$\text{hi}(k) = \min\bigl(\text{hi}(k-1),\; \Sigma_k + \text{uncertainty}(k)\bigr)$$

where $\Sigma_k = \sum_{j=1}^{k} \Delta E_j$ is the cumulative energy sum after step $k$. The ratchet ensures lo is non-decreasing and hi is non-increasing across steps. Because $\text{uncertainty}(k) = \varepsilon_0/\sqrt{k}$ is strictly decreasing, the *proposed* interval shrinks each step and the max/min operators guarantee the stored interval never widens.

The initial conditions depend on the constructor used:
- **Without pre-seed**: step 1 sets $[\Sigma_1 - \varepsilon_0, \Sigma_1 + \varepsilon_0]$ directly
- **With pre-seed** (Phase 6 constructor): interval starts at $[\text{MJ}_\text{direct} - \varepsilon_0, \text{MJ}_\text{direct} + \varepsilon_0]$; subsequent steps ratchet inward from there

Crossover ($\text{lo}(k) > \text{hi}(k)$) is detected and the interval resets to the unconstrained proposed bounds for that step, flagging `CrossoverOccurred` for diagnostic logging.

### 3.3 Contact Shell Cutoff

BFS terminates when all residues within the **8 Å Cα–Cα shell** around $v_m$ have been visited (Phase 1). In practice, the BFS visits every connected residue in the graph reachable from $v_m$ by the contact edge types in §1.2, stopping when no new residues contribute energy changes above a numerical floor ($|\Delta E_k| < 10^{-9}$ kcal/mol).

---

## 4. ε₀ Calibration

$\varepsilon_0$ is the single free parameter of the MutCert model. It controls the initial interval half-width and is calibrated once per protein family from the training split.

### 4.1 Calibration Procedure

1. Run the BFS agent on all training mutations $\mathcal{T}$ (80% split, stratified by position)
2. Compute the OLS linear recalibration: $\widehat{\Delta\Delta G} = \alpha + \beta \cdot \text{MJ}_{\text{direct}}$
3. Compute residuals: $r_i = \Delta\Delta G_i^{\exp} - \widehat{\Delta\Delta G}_i$
4. Set $\varepsilon_0 = \text{RMSE}(\mathcal{T}) = \sqrt{\frac{1}{|\mathcal{T}|} \sum_{i} r_i^2}$

This makes the interval half-width at hop 0 equal to the training RMSE — a principled, data-driven choice that avoids arbitrary thresholds.

### 4.2 Observed Values

| Family | Train n | RMSE (train) | $\varepsilon_0$ | Spearman ρ (val) |
|--------|:-------:|:------------:|:---------------:|:----------------:|
| T4 lysozyme | 97 | ~1.79 kcal/mol | 1.7886 | −0.449 |
| CI2 | 27 | ~0.73 kcal/mol | (see run output) | −0.194 |
| Barnase | 46 | ~1.12 kcal/mol | (see run output) | −0.441 |

*Exact ε₀ values are printed by `validate-all` and stored in run-registry per family.*

---

## 5. Linear Recalibration (OLS)

After the BFS agent computes the raw MJ direct energy sum $\text{MJ}_{\text{direct}}$ for each mutation, a linear recalibration corrects for systematic bias:

$$\widehat{\Delta\Delta G} = \hat{\alpha} + \hat{\beta} \cdot \text{MJ}_{\text{direct}}$$

where $\hat{\alpha}$ and $\hat{\beta}$ are estimated by ordinary least squares on the training set. This corrects:
- Global scale bias (MJ units vs kcal/mol)
- Additive offset (WT stability baseline)

The recalibrated estimate is used as the interval midpoint. The interval width is determined by the max/min ratchet (§3.2) independently of the OLS fit.

**Implementation**: `ValidationPipeline.cs` — step 3.7 (OLS fit on `trainGated`), step 4 (val loop with recalibrated DDG).

---

## 6. Wild-Type Gate

Before OLS fitting and validation, a **WT gate** filters records where the graph residue type does not match the wild-type code in S2648. This can arise from:
- Residual sequence offset errors (e.g., CIF propeptide residues at wrong position)
- AlphaFold prediction artefacts at termini
- S2648 numbering inconsistencies

Filtered records are counted and logged but excluded from both training and validation. This ensures no energy calculations are performed on "wrong" residues that would inject noise.

**Effect (Phase 9-E)**: Barnase Spearman ρ improved from −0.397 to −0.441 after gating 3 mismatched validation records.

**Implementation**: `ValidationPipeline.cs` — step 3.5 (WT audit), step 3.6 (WT gate producing `trainGated` / `valGated`).

---

## 7. Summary of Hyperparameters

| Parameter | Value | Set by |
|-----------|:-----:|--------|
| Contact shell radius | 8 Å (VdW/hydrophobic) or 10 Å (electrostatic) | Code (StructureLoader.cs) |
| Interval rule | max/min ratchet, uncertainty = ε₀/√k | Code (EnergySignalPropagator.cs) |
| Width threshold $\theta_w$ | 2.0 kcal/mol | Design (Phase 1) |
| Monotone window $\ell$ | 5 steps | Design |
| Minimum steps $k_{\min}$ | 5 | Design |
| Train/val split | 80/20 | Fixed (stratified) |
| $\varepsilon_0$ | Per-family RMSE | Calibrated |
| OLS intercept $\hat{\alpha}$ | Per-family | Calibrated |
| OLS slope $\hat{\beta}$ | Per-family | Calibrated |

**Note**: MutCert has **zero** neural network parameters and **zero** tunable force-field weights. The only calibrated parameters are $\varepsilon_0$, $\hat{\alpha}$, $\hat{\beta}$ — all estimated in closed form from the training set.

---

## 8. Reproducibility

All results are fully reproducible:
1. Frozen graphs in Neo4j (immutable after import)
2. Deterministic BFS order (graph creation order)
3. Deterministic OLS (closed-form solution)
4. All RunStep records stored in `run-registry` (glass-box trace available for every run)
5. Seed-free — no stochastic elements

To reproduce Phase 9-E results:
```
dotnet run --project src/Protein.Engine -- validate-all
```

---

*Generated: Phase 11 — Preprint Prep*  
*Keywords: methods-section, miyazawa-jernigan-1996, reference-state, bfs-propagation, ratchet-narrowing, epsilon0-calibration, ols-recalibration, wt-gate, reproducibility, phase11*
