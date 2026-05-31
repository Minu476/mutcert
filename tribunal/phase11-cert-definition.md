# MutCert — Convergence Certificate: Formal Definition

**Phase 11 Preprint Prep | MutCert v1.0**  
**Keywords:** convergence-certificate, formal-definition, latex, monotone-narrowing, width-threshold, ratchet-rule, phase11

---

## 1. Informal Statement

A MutCert Convergence Certificate is a **machine-issued record** that the energy-interval estimate $[\text{lo}(k), \text{hi}(k)]$ for a given point mutation has narrowed to a pre-specified precision threshold and has done so monotonically.

The certificate records *precision*, not accuracy: the interval is guaranteed to be narrow, but is not guaranteed to contain the true ΔΔG. The accuracy of the midpoint is assessed separately via Spearman ρ on the validation set. See §6 for the known limitation illustrated by the L121A trace.

The certificate is stored as a `ConvergenceCertificate` node in the `run-registry` Neo4j database, linked to every BFS step that contributed to it.

---

## 2. Definitions

| Symbol | Type | Definition |
|--------|------|------------|
| $m$ | String | Mutation identifier, e.g., `MUT_P00720_A_121_ALA` |
| $G = (V, E)$ | Graph | Frozen residue contact graph for the protein family |
| $v_m \in V$ | Node | Mutation site residue node |
| $k$ | Integer ≥ 1 | BFS step index (1-indexed) |
| $\text{lo}(k), \text{hi}(k)$ | $\mathbb{R}$ | Lower/upper bounds of the ΔΔG interval at step $k$ |
| $w(k)$ | $\mathbb{R}_{\geq 0}$ | Interval width: $w(k) = \text{hi}(k) - \text{lo}(k)$ |
| $\Delta E_k$ | $\mathbb{R}$ | Energy contribution of step $k$ (MJ contact energy × hop decay) |
| $\varepsilon_0$ | $\mathbb{R}_{> 0}$ | Calibrated interval half-width parameter |
| $\theta_w$ | $\mathbb{R}_{> 0}$ | Width threshold for certificate issuance |
| $\ell$ | Integer ≥ 1 | Monotone window length (consecutive non-increasing steps required) |
| $k_{\min}$ | Integer ≥ 1 | Minimum BFS steps before certificate can be issued |

---

## 3. Formal Certificate Conditions

### 3.1 LaTeX Source

```latex
\begin{definition}[Convergence Certificate]
Let $m$ be a point mutation on protein graph $G = (V, E)$,
and let $\{[\mathrm{lo}(k), \mathrm{hi}(k)]\}_{k=1}^{K}$ be the
sequence of interval estimates produced by the FSDE BFS agent at each
step $k$. Define interval width $w(k) := \mathrm{hi}(k) - \mathrm{lo}(k)$.

A \emph{Convergence Certificate} $\mathcal{C}(m, k^*)$ is issued at the
first step $k^* \geq k_{\min}$ such that the following three conditions
hold simultaneously:

\begin{align}
  \text{(C1)} \quad & w(k^*) \leq \theta_w
    & \text{[precision threshold]} \\[4pt]
  \text{(C2)} \quad & w(k^* - j + 1) \leq w(k^* - j)
    \;\; \forall j \in \{1, \ldots, \ell - 1\}
    & \text{[monotone narrowing over } \ell \text{ steps]} \\[4pt]
  \text{(C3)} \quad & k^* \geq k_{\min}
    & \text{[minimum shell coverage]}
\end{align}

The certificate $\mathcal{C}(m, k^*)$ records:
\begin{itemize}
  \item the mutation identifier $m$,
  \item the certified interval $[\mathrm{lo}(k^*), \mathrm{hi}(k^*)]$,
  \item the midpoint $\widehat{\Delta\Delta G} = \tfrac{\mathrm{lo}(k^*) + \mathrm{hi}(k^*)}{2}$,
  \item the step index $k^*$ at which the certificate was issued,
  \item and the full ordered sequence of RunStep records
        $(v_k, \Delta E_k, \mathrm{hop}_k, \mathrm{lo}(k), \mathrm{hi}(k))_{k=1}^{k^*}$
        (the \emph{glass-box causal chain}).
\end{itemize}
\end{definition}
```

### 3.2 Parameter Values (Phase 1, v1.0)

| Parameter | Symbol | Phase 1 value | Phase 2 target |
|-----------|--------|:-------------:|:--------------:|
| Width threshold | $\theta_w$ | 2.0 kcal/mol | 1.0 kcal/mol |
| Monotone window | $\ell$ | 5 consecutive steps | 5 consecutive steps |
| Minimum steps | $k_{\min}$ | 5 | 8 |
| Contact shell radius | $r_{\text{shell}}$ | 8 Å | 12 Å |

---

## 4. Interval Update Rule

The interval is updated by the **max/min ratchet** implemented in `EnergySignalPropagator.cs`:

$$\text{uncertainty}(k) = \varepsilon_0 / \sqrt{k}$$
$$\text{proposed\_lo}(k) = \Sigma_k - \text{uncertainty}(k)$$
$$\text{proposed\_hi}(k) = \Sigma_k + \text{uncertainty}(k)$$
$$\text{lo}(k) = \max\bigl(\text{lo}(k-1),\; \text{proposed\_lo}(k)\bigr)$$
$$\text{hi}(k) = \min\bigl(\text{hi}(k-1),\; \text{proposed\_hi}(k)\bigr)$$

where $\Sigma_k = \sum_{j=1}^{k} \Delta E_j$ is the cumulative energy sum after $k$ BFS steps.

**Why this rule satisfies C2 by construction**: because $\varepsilon_0/\sqrt{k}$ is monotone decreasing, each proposed interval is strictly narrower than the previous one. The max/min operators then guarantee the stored bounds never move outward. Monotone narrowing is therefore a structural consequence of the algorithm, not an empirical property that can fail.

**What C2 does not guarantee**: the midpoint $\frac{\text{lo}(k)+\text{hi}(k)}{2}$ may be far from the true ΔΔG even when the interval is narrow. The interval converges to the running cumulative MJ sum, which is only a good estimate when the MJ potential is accurate for that mutation.

Initial conditions: without pre-seed, step 1 sets $[\Sigma_1 - \varepsilon_0, \Sigma_1 + \varepsilon_0]$ directly. With pre-seed, the interval is initialised at $[\text{MJ}_\text{direct} - \varepsilon_0, \text{MJ}_\text{direct} + \varepsilon_0]$.

Crossover recovery: if $\text{lo}(k) > \text{hi}(k)$, the interval resets to the unconstrained proposed bounds and `CrossoverOccurred` is flagged.

---

## 5. Calibration of $\varepsilon_0$

$\varepsilon_0$ is the primary calibration parameter. It is computed once per protein family from the held-out training set:

$$\varepsilon_0 = \sqrt{\frac{1}{|\mathcal{T}|} \sum_{i \in \mathcal{T}} \left(\Delta\Delta G_i^{\exp} - \widehat{\Delta\Delta G}_i\right)^2}$$

where $\mathcal{T}$ is the training split and $\widehat{\Delta\Delta G}_i$ is the midpoint of the MutCert interval for mutation $i$ **before** applying the linear recalibration offset.

This makes the initial interval width $2\varepsilon_0$ equal to the root-mean-square residual error on the training set — a principled, data-driven choice.

**Observed values (Phase 9-E run)**:
- T4 lysozyme: $\varepsilon_0 \approx 1.789$ kcal/mol
- CI2: $\varepsilon_0 \approx 0.xx$ kcal/mol (see run output)
- Barnase: $\varepsilon_0 \approx 0.xx$ kcal/mol (see run output)

---

## 6. Known Limitation: Precision Without Accuracy

The certificate certifies that the interval is narrow and narrowed monotonically. It does not certify that the interval contains the true ΔΔG.

**Illustrative case — MUT_P00720_A_121_ALA** (T4 lysozyme, Phase 9-E run):
- Certificate ISSUED at step 53
- Final interval: [−2.47, −1.98] kcal/mol, width = 0.49 kcal/mol (≤ 2.0, C1 satisfied)
- Experimental ΔΔG: approximately +22.8 kcal/mol (destabilising)
- The certified interval is ~25 kcal/mol away from the experimental value

This happens because the MJ potential assigns negative (stabilising) contact-energy change for alanine substitutions at buried hydrophobic positions, while the real effect is large destabilisation. The ratchet faithfully narrows around the running MJ sum; the MJ sum is wrong for this mutation.

The correct framing for a preprint is: MutCert issues a **convergence certificate for the interval produced by the MJ-BFS estimate**. Whether that estimate is accurate depends on the quality of the underlying energy function, which is a separate question assessed by Spearman ρ on the validation set.

A certificate that is systematically miscalibrated for large destabilising buried mutations is a concrete research finding, not a failure to hide.

---

## 7. What the Certificate Is and Is Not

| Property | MutCert v1 |
|----------|:-----------:|
| Interval narrows monotonically | **by construction** (ratchet with ε₀/√k) |
| Certificate conditions inspectable | ✓ |
| Full BFS causal chain stored and auditable | ✓ |
| Certified interval contains true ΔΔG | **not guaranteed** |
| Accuracy assessed by | Spearman ρ on held-out val set |

DDGun, FoldX, and Rosetta all produce point estimates without convergence tracking. MutCert adds the tracking layer. The value of the tracking layer is that it makes the *convergence behaviour* of the energy-sum algorithm fully visible, which is useful for diagnosing where the MJ estimate goes wrong.

---

*Generated: Phase 11 — Preprint Prep*  
*Keywords: convergence-certificate, formal-definition, c1-width-threshold, c2-monotone-narrowing, c3-minimum-steps, epsilon0-calibration, ratchet-rule, precision-not-accuracy, latex, phase11*
