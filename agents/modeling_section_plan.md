# Modelling Section Plan

> Section 2 of the JDSMC paper: **"The SCR-ASC System: Mathematical Modelling"**

This document plans the content, structure, and sourcing for the Modelling section. The goal is to condense the thesis's Chapter 3 (Governing Equations) and parts of Chapter 4 (NOx Model) into a self-contained journal section that derives the discrete nonlinear switched model from first principles.

---

## Source Material (Thesis)

| Thesis Source | Content | Journal Target |
|---|---|---|
| `Thesis/secs/3-governing_eqns/1-reactions.tex` | Full reaction list, reduced 3 reactions, Arrhenius/polynomial rate constant model, ASC lumping | §2.1, §2.2 |
| `Thesis/secs/3-governing_eqns/2-NH3ads_proc_dyn.tex` | NH₃ adsorption/desorption rate eqns, molar conservation at residence-time scale, telescopic summation, σ(k) derivation | §2.3 |
| `Thesis/secs/3-governing_eqns/2-1-sig_approx.tex` | Zero-order hold approximation of σ(k), introduction of g_sat switching function | §2.3, §2.5 |
| `Thesis/secs/3-governing_eqns/3-NOx_proc_dyn.tex` | NOₓ process dynamics, molar conservation, ZOH + avg surface conc approximation | §2.4 |
| `Thesis/secs/3-governing_eqns/4-NH3_proc_dyn.tex` | Gaseous NH₃ dynamics (note: excluded from journal I/O model per thesis §6) | Omit or brief mention |
| `Thesis/secs/3-governing_eqns/5-urea_dosing_proc_dyn.tex` | Urea dosing → inlet NH₃ model | §2.2 |
| `Thesis/secs/3-governing_eqns/6-proc_dyn.tex` | Summary of complete process dynamics (4 eqns) | §2.5 |
| `Thesis/secs/4-NOx_mdl/1_NOx_r_mdl.tex` | σ elimination → η(k) recursive model, parametrization, identifiable form | §2.6 |
| `Thesis/secs/4-NOx_mdl/2_catalyst_saturation.tex` | Saturated catalyst model, min-max form, LP for θ_sat | §2.7 |
| `Thesis/secs/4-NOx_mdl/25_switching_condition.tex` | Switching condition parametrization, guard condition, SETAR structure | §2.6–§2.7 |
| `Thesis/secs/4-NOx_mdl/1-subs/parametrization.tex` | η dynamics parametrization, Chebyshev basis selection rationale (φ₁ vs φ₂), assumptions A3–A6 | §2.2, §2.6 |
| `Thesis/secs/2-data/1-res_time.tex` | Residence time definition, RTD mode, τ_mode ≈ 0.07 s *(chamber dimensions → validation section)* | §2.2 |
| `Thesis/secs/2-data/2-data_units.tex` | Exhaust gas density (ideal gas law), constant-density justification, density sensitivity to NOₓ *(ppm ↔ mol/m³ conversion → validation section)* | §2.2 |

---

## Proposed Subsection Structure

### §2.1 — SCR-ASC Reactions and Modelling Assumptions
**File:** `secs/2-modeling/1-reactions.tex`

- List the full set of reactions occurring in the SCR-ASC chamber (from thesis §3.1, abridged).
- State the **three dominant reactions** retained for reduced-order modelling:
  1. Standard SCR: $4 NH_3^{ads} + 4 NO + O_2 \xrightarrow{k_{scr}} 4 N_2 + 6 H_2O$
  2. Ammonia Oxidation (AMOX): $4 NH_3^{ads} + 3 O_2 \xrightarrow{k_{oxi}} 2 N_2 + 6 H_2O$
  3. Adsorption/Desorption: $NH_3 + \Theta_{free} \xrightleftharpoons[k_{des}]{k_{ads}} NH_3^{ads}$
- **Assumptions** (consolidated list):
  - Eley-Rideal mechanism for SCR reaction
  - All NOₓ is NO (sensor indistinguishable)
  - Slow SCR reaction omitted (flow-rate dominated)
  - Mass transfer neglected (reaction-controlled)
  - 100% nitrogen selectivity for AMOX (both SCR and ASC)
  - ASC AMOX lumped into SCR AMOX (single oxidation reaction with combined rate constant)
- Briefly state the reaction rate expressions and units (surface vs. volumetric).
- **Figure:** Schematic of the reduced SCR-ASC model (thesis Fig 3.2).

> **Note:** The intro section already presents the CSTR model and its limitations. This subsection should **not** repeat the CSTR discussion but rather start directly with the discrete modelling approach.

---

### §2.2 — Physical Property Models
**File:** `secs/2-modeling/2-physical_models.tex`

This subsection collects all the approximate models for physical properties used to parametrize the governing equations. These models bridge the gap between the first-principles rate equations (§2.1) and the identifiable parametric form (§2.6).

#### 2.2.1 — Exhaust Gas Density
*(Source: `Thesis/secs/2-data/2-data_units.tex`, §2.1)*

- Exhaust gas density approximated via **ideal gas law** (assuming air-like composition at ambient pressure):
  $$\rho = \frac{PM}{RT} = \frac{\mu}{T}$$
  where $P = 101.325$ kPa, $M = 28.97$ g/mol, $R = 8.314$ J/(mol·K).
- **Constant-density justification (Assumption A6):** At the typical operating temperature ($T_0 \approx 250°C \approx 523$ K), a ±100°C variation yields < 10% change in density — negligible compared to mass flow rate variations. Thus:
  $$\rho \approx \rho_0 = 1.2 \times 10^{-3} \, \text{g/cm}^3 \quad \text{(at } 250°C\text{)}$$
- **Density sensitivity to NOₓ concentration:** Also shown to be negligible because the mole fraction of NOₓ in exhaust is extremely small (ppm-level).

#### 2.2.2 — Residence Time Model
*(Source: `Thesis/secs/2-data/1-res_time.tex`, `Thesis/secs/4-NOx_mdl/1-subs/parametrization.tex` — Assumption A6)*

- **Definition:** $\tau = \rho V_{scr} / F$ — the time a fluid parcel remains in the control volume.
- Using constant density: $\tau(k) = \rho_0 V_{scr} / F(k) = \tau_0 / F(k)$.
- **Residence time distribution (RTD):** Due to large flow-rate variations, the **mode** of the RTD is used:
  $$\tau_{mode}^{test} \approx 0.072 \, \text{s}, \quad \tau_{mode}^{truck} \approx 0.071 \, \text{s}$$
- This is much shorter than sampling intervals ($t_s^{test} = 0.2$ s, $t_s^{truck} = 1$ s), confirming $n = t_s/\tau \gg 1$ — the fundamental motivation for the telescopic summation approach.

> *Chamber dimensions table and ppm ↔ mol/m³ unit conversions deferred to the [validation section plan](file:///mnt/c/Users/sesha/My%20Drive/scharla/ResearchWork/SCR-OBD/Publications/SCR-dynamics_JDSMC/agents/validation_section_plan.md).*

#### 2.2.3 — Rate Constant Temperature Models (Arrhenius → Polynomial)
*(Source: `Thesis/secs/3-governing_eqns/1-reactions.tex`, §3.1.1–3.1.2)*

- Arrhenius equation: $k = A \exp(-E/RT)$.
- **Linear approximation** (1st-order Taylor): $k(T) \approx mT + c$, valid within ±50°C of $T_0$.
- **Quadratic approximation** (2nd-order Taylor): $k(T) \approx qT^2 + mT + c$, for wider temperature range.
- **Chebyshev basis** (from `parametrization.tex`): For numerical stability, the temperature range $[T_{min}, T_{max}]$ is mapped to $[-1, 1]$:
  $$\phi_1(k) = \begin{bmatrix} \frac{T-T_0}{T_r} & 1 \end{bmatrix}, \quad \phi_2(k) = \begin{bmatrix} 2\left(\frac{T-T_0}{T_r}\right)^2 - 1 & \frac{T-T_0}{T_r} & 1 \end{bmatrix}$$
  where $T_0 = (T_{max}+T_{min})/2$ and $T_r = (T_{max}-T_{min})/2$.
- **Basis order selection rationale** (from `parametrization.tex`, assumptions A3–A4):
  - $\phi_1$ (linear) is used for terms containing **individual rate constants** ($k_{ads}$, $k_{od}$, $k_{scr}$), which vary nearly linearly over the operating range.
  - $\phi_2$ (quadratic) is used for terms containing **products of rate constants and Γ** (e.g., $k_{scr} \cdot k_{ads} \cdot \Gamma$), where the combined temperature dependence includes an inflection point.

#### 2.2.4 — Urea Dosing → Inlet NH₃ Model
*(Source: `Thesis/secs/3-governing_eqns/5-urea_dosing_proc_dyn.tex` — Assumption A5)*

- Urea (AdBlue, 32.5% aqueous urea solution) decomposes to ammonia. The aggregated reaction:
  $NH_2\text{-}CO\text{-}NH_2 \longrightarrow 2 NH_3 + CO_2 + x H_2O$
- Rate is constant (fixed urea concentration in solution). Inlet NH₃ concentration:
  $$[NH_3]^{in}(k) = \nu_u \frac{u_{inj}(k)}{F(k)}$$
  where $\nu_u = 2 r_u t_s \rho_0$ is a lumped constant (Assumption A5: independent of temperature since urea is preheated).
- The reciprocal dependence on $F$ breaks down at zero flow (physically consistent — no flow means no reaction).

> **Decision needed:** Should we use linear or quadratic Chebyshev basis for the rate constants in this paper? The thesis uses both; the parametrization file clarifies the rationale (linear for single rate constants, quadratic for rate–Γ products). This rationale should be stated explicitly in §2.2.5.

---

### §2.3 — Ammonia Adsorption/Desorption Process Dynamics
**File:** `secs/2-modeling/3-NH3ads_dynamics.tex`

- **Rate equations** for adsorption/desorption (surface concentrations):
  $$\dot{[NH_3]}^{ads} = \Gamma \gamma_{ads}(k) - [NH_3]^{ads} \gamma_{des}(k)$$
  with $\gamma_{ads} = k_{ads}[NH_3]^{in}$ and $\gamma_{des} = k_{ads}[NH_3]^{in} + k_{des} + k_{scr}[NO_x]^{in} + k_{oxi}$.
- **Molar conservation at residence-time scale** (thesis §3.2.1):
  - Within $t_s$, there are $n = t_s/\tau$ residence times.
  - Telescopic summation yields $\Omega(k)$, the total surface concentration change.
  - Define $\sigma(k) = \frac{1}{n}\sum_{i=0}^{n-1} [NH_3]^{ads}(k + i\tau)$ — the **average surface concentration**.
- **Zero-order hold approximation:** $\sigma(k) \approx [NH_3]^{ads}(k)$ leading to the discrete update:
  $$\sigma^{ub}(k+1) = \sigma(k) + t_s k_{ads}[NH_3]^{in}(\Gamma - \sigma(k)) - t_s(k_{oxi} + k_{des})\sigma(k) - t_s k_{scr}[NO_x]^{in}\sigma(k)$$
- **Figure:** Discrete plug-flow reactor model (thesis Fig 3.1).

---

### §2.4 — NOₓ Process Dynamics
**File:** `secs/2-modeling/4-NOx_dynamics.tex`

- **Rate equation:** $\dot{[NO_x]}^{scr} = -k_{s2v} k_{scr} [NH_3]^{ads} [NO_x]^{in}$
- **Molar conservation** at residence-time scale (note: **no integrating effect** — outlet depends on inlet one τ earlier).
- **ZOH + average surface conc approximation** yields:
  $$[NO_x]^{out}(k+1) = [NO_x]^{in}(k)(1 - k_{s2v} k_{scr} \sigma(k) \tau)$$
- This is a key result: the NOₓ output is a **static map** modulated by σ(k).

---

### §2.5 — Complete Process Dynamics and Catalyst Saturation
**File:** `secs/2-modeling/5-process_dynamics.tex`

- Collect the four governing equations (thesis §3.6):
  1. Urea dosing: $[NH_3]^{in}(k)$
  2. NOₓ reduction: $[NO_x]^{out}(k+1)$
  3. Gaseous NH₃ (mention for completeness; not used for I/O model)
  4. Ammonia adsorption/desorption: $\sigma^{ub}(k+1)$ with $g_{sat}$
- Introduce the **switching function** $g_{sat}(\sigma)$ (already defined in intro §1, eqn. (2)).
- State: tailpipe NH₃ dynamics don't affect NOₓ reduction and tailpipe ammonia measurements are unavailable → focus on NOₓ I/O model.

> **Note:** This subsection is a brief transition paragraph collecting the process dynamics before the σ-elimination step. It can be kept concise — possibly merged with §2.4 or §2.6.

---

### §2.6 — NOₓ Reduction Input-Output Model: σ Elimination and Parametrization
**File:** `secs/2-modeling/6-NOx_IO_model.tex`

This is the **core contribution** of the modelling section.

- **Define** $\eta(k) = u_1(k-1) - x_1(k)$ (fractional NOₓ reduction).
- **σ-elimination:** Using the NOₓ output equation $\sigma(k) = \eta(k+1) / (\tau k_{s2v} k_{scr} u_1(k))$, substitute into the σ update equation to obtain the **recursive η model** (thesis eqn. (4.8)):
  $$\eta(k+1) = \eta(k) \cdot f_{flow}(k) \cdot [1 - t_s(k_{ads}[NH_3]^{in} + k_{scr}u_1 + k_{od})] + t_s k_{s2v} \Gamma \tau u_1 [NH_3]^{in} k_{scr} k_{ads}$$
- **Parametrization** with Chebyshev temperature models:
  - Introduce $\bar{\eta}_F(k) = F(k-1)\eta(k)/u_1(k-1)$ (flow-scaled fractional NOₓ reduction).
  - Four parameter groups: $\theta_\Gamma$, $\theta_{\eta_{ads}}$, $\theta_{\eta_{od}}$, $\theta_{\eta_{scr}}$.
  - Compact linear-in-parameters form: $\bar{\eta}_F(k+1) = \bar{\eta}_F(k) + \phi_{NO_x}^T \theta_{NO_x}$.
- **Physical interpretation** of each parameter group.
- Discuss the **identifiable form** (thesis §4.1.3): the regressor $\phi_{NO_x}$ depends only on measured quantities.

---

### §2.7 — Switched Model: Saturation, Guard Conditions, and Bounds
**File:** `secs/2-modeling/7-switched_model.tex`

- **Saturated catalyst dynamics:** When $\sigma = \Gamma$, adsorption dynamics switch off → derive $f_\Gamma(k) = (u_1(k)/F(k)) \phi_2^T(k) \theta_{sat}$.
- **Guard condition** (thesis eqn. (4.25)): The linear-in-parameters switching condition $\phi_g^T \theta_{NO_x} \leq 0$ determines the active mode.
- **Min-max formulation:** $\eta(k+1) = \max\{0, \min\{f_\sigma(k), f_\Gamma(k)\}\}$.
- **Bounds on tailpipe NOₓ:** $\max\{0, u_1(k) - f_\Gamma(k)\} \leq x_1(k+1) \leq u_1(k)$.
- **LP for θ_sat estimation:** Minimize saturated-mode NOₓ reduction subject to $f_\Gamma(k) \geq \eta(k+1)$ for all $k$.
- **State diagram** figure showing the switching between saturated/unsaturated modes.
- Identify the model as a **self-excited threshold nonlinear autoregressive with exogenous input (SETAR)** structure.

---

## Figures Needed

| # | Description | Source | Status |
|---|---|---|---|
| 1 | SCR-ASC reduced model schematic (lumped ASC into SCR) | Thesis Fig 3.2 | ☐ Adapt for journal |
| 2 | Discrete plug-flow reactor model | Thesis Fig 3.1 | ☐ Adapt for journal |
| 3 | SCR-ASC system abstraction (control volume) | Thesis Fig 4.1 | ☐ Adapt for journal |
| 4 | State diagram for NOₓ switching model | Thesis Fig 4.4 | ☐ Adapt for journal |

> **Tip:** Figures 2 and 3 may be combinable into a single figure for the journal.

---

## LaTeX File Plan

| File | Content |
|---|---|
| `secs/2-modeling/0-modeling.tex` | Section header + brief preamble + `\input` calls |
| `secs/2-modeling/1-reactions.tex` | §2.1 — Reactions and assumptions |
| `secs/2-modeling/2-physical_models.tex` | §2.2 — Rate constant models, residence time, urea dosing |
| `secs/2-modeling/3-NH3ads_dynamics.tex` | §2.3 — Adsorption/desorption dynamics |
| `secs/2-modeling/4-NOx_dynamics.tex` | §2.4 — NOₓ process dynamics |
| `secs/2-modeling/5-process_dynamics.tex` | §2.5 — Summary + saturation |
| `secs/2-modeling/6-NOx_IO_model.tex` | §2.6 — I/O model derivation (core contribution) |
| `secs/2-modeling/7-switched_model.tex` | §2.7 — Switched model and estimation |

---

## Open Questions

1. **Gaseous NH₃ dynamics:** The thesis explicitly excludes these from the I/O model (§6, line 44–46). Should we omit §2.4-style derivation for NH₃ entirely, or include a brief remark for completeness?

2. **Chebyshev basis order:** Linear ($\phi_1$) for rate constants and quadratic ($\phi_2$) for rate–Γ products? Or simplify to just linear throughout?

3. **Level of detail for σ derivation:** The telescopic summation (thesis §3.2.1–3.2.3) is about 2.5 pages. For the journal, should we present the full derivation or summarize the key steps and point to a reference/appendix?

4. **Validation figure from thesis Fig 4.5 (3_validation.tex):** Should model validation be placed here in §2 or deferred to the System Identification section (§3)?
