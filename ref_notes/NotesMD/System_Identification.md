# System Identification for the NOx Dynamic Model

> Continuation of [NOx_Dynamic_Model_Derivation.md](NOx_Dynamic_Model_Derivation.md).  
> Equations (S), (N), (U), (E), (SW), (D) and notation ($\phi_1$, $\phi_2$, $\phi_{NO_x}$, $\theta_{NO_x}$, $\theta_{sat}$, $G(k)$, $\bar{\eta}_F$) are defined there.  
> Source: Thesis Chapters 2, 4, 5.

---

The switched model (Section 4d of the derivation) has two parameter vectors: $\theta_{sat}$ (saturation ceiling, Section 4b) and $\theta_{NO_x}$ (unsaturated dynamics and guard condition, Sections 4a/4c). The goal is to estimate both from available measurement data ($u_1, x_1, u_2, F, T$). The guard condition $G(k) = \phi_g^T\,\theta_{NO_x}$ depends on $\theta_{NO_x}$, so the operating mode is not known a priori. We address this by first estimating $\theta_{sat}$ via LP (mode-independent), then using it to identify saturated segments, and finally estimating $\theta_{NO_x}$ from the remaining data.

---

## 1. Saturated model parameters ($\theta_{sat}$)

### 1a. Bounding condition

From Section 4b, $\bar{\eta}_{F,sat}(k) = \phi_2^T(k)\,\theta_{sat}$. This ceiling bounds the actual reduction at every time step:

$$
\eta_{sat}(k) = \phi_{sat}^T(k)\,\theta_{sat} \ge \eta(k) \qquad \forall\, k
$$

where $\phi_{sat}(k) = \frac{u_1(k)}{F(k)}\,\phi_2(T(k))$ is the saturated-model regressor (converting $\bar{\eta}_{F,sat}$ from the derivation back to $\eta$ space). This holds because $\sigma = \Gamma$ always produces at least as much reduction as any $\sigma \le \Gamma$. The bound is mode-independent, so all data can be used.

### 1b. LP formulation

Since the bound must hold for every time step, the optimal $\hat{\theta}_{sat}$ is the one that produces the tightest possible envelope over the data. This is a linear program: minimize the total area under the saturated curve subject to the bound constraint:

$$
\hat{\theta}_{sat} = \arg\min_{\theta_{sat}} \;\mathbf{1}^T\,\Phi_{sat}\,\theta_{sat} \qquad \text{s.t.} \quad \Phi_{sat}\,\theta_{sat} \succeq H
$$

where $\Phi_{sat} = [\phi_{sat}(1), \ldots, \phi_{sat}(N{-}1)]^T$ and $H = [\eta(2), \ldots, \eta(N)]^T$.

The operating temperature range is partitioned into two zones to capture the non-monotone dependence of $k_{scr}\,\Gamma$ on $T$. "High" refers to the upper partition (RMC: 250-350°C; hot-FTP: 200-300°C); "low" refers to the lower partition (cold-FTP only, below 200°C). Each zone has its own $\theta_{sat}$.

| Test | Temp. zone | DG $\theta_{sat}[2]$ | DG $\theta_{sat}[1]$ | DG $\theta_{sat}[0]$ | Aged $\theta_{sat}[2]$ | Aged $\theta_{sat}[1]$ | Aged $\theta_{sat}[0]$ |
|------|-----------|---|---|---|---|---|---|
| RMC | high | -0.06 | 1.43 | 31.19 | -0.07 | 1.77 | 27.82 |
| hot-FTP | high | -0.49 | 2.94 | 40.94 | -0.54 | 3.40 | 39.59 |
| cold-FTP | high | -0.19 | 0.97 | 40.69 | -0.28 | 1.82 | 38.97 |
| cold-FTP | low | 0.26 | 5.97 | 45.00 | 0.39 | 8.25 | 50.63 |

Bounded $\eta$ plots (saturated ceiling as tight envelope over actual data):

| Degreened | Aged |
|-----------|------|
| ![DG cold-FTP](../Thesis/figs/4-NOx_mdl/2_figs/bounded_eta_plots/eta_bounds_dg_cftp.png) | ![Aged cold-FTP](../Thesis/figs/4-NOx_mdl/2_figs/bounded_eta_plots/eta_bounds_aged_cftp.png) |
| ![DG hot-FTP](../Thesis/figs/4-NOx_mdl/2_figs/bounded_eta_plots/eta_bounds_dg_hftp.png) | ![Aged hot-FTP](../Thesis/figs/4-NOx_mdl/2_figs/bounded_eta_plots/eta_bounds_aged_hftp.png) |
| ![DG RMC](../Thesis/figs/4-NOx_mdl/2_figs/bounded_eta_plots/eta_bounds_dg_rmc.png) | ![Aged RMC](../Thesis/figs/4-NOx_mdl/2_figs/bounded_eta_plots/eta_bounds_aged_rmc.png) |

---

## 2. Segment identification using $\theta_{sat}$

With $\hat{\theta}_{sat}$ estimated, the saturated response $\hat{\eta}_{sat}(k) = \phi_{sat}^T(k)\,\hat{\theta}_{sat}$ is computed for each time step. A data point is classified as near-saturated if:

$$
|\hat{\eta}_{sat}(k) - \eta(k)| \le \varepsilon_{sat} \implies \text{saturated operation}
$$

The threshold $\varepsilon_{sat}$ is chosen from the prediction error variance ($\varepsilon_{sat} = 2.5 \times 10^{-3}\;\text{mol/m}^3$ for the test-cell dataset). The saturated response depends only on $u_1$, $F$, $T$ (not on the previous state or urea dosing), so it can be evaluated without simulation.

![Catalyst mode detection for DG RMC](../Thesis/figs/5-sat_detect/3_parm_ID/SatDetect_dg_rmc_1_FTIR.png)

| Test | $N$ | $N_{sat}$ (FTIR) | $N_{sat}$ (NOx sensor) |
|------|-----|-------------------|------------------------|
| DG RMC 1 | 2402 | 490 | 492 |
| DG RMC 2 | 2402 | 493 | 494 |
| DG RMC 3 | 2402 | 491 | 493 |
| Aged RMC | 2402 | 515 | 515 |
| DG Hot FTP 1 | 634 | 186 | 194 |
| DG Hot FTP 2 | 634 | 181 | 191 |
| DG Hot FTP 3 | 638 | 180 | 198 |
| Aged Hot FTP | 636 | 198 | 220 |

The unsaturated segments (the complement) are used for $\theta_{NO_x}$ estimation next.

---

## 3. Unsaturated model parameters ($\theta_{NO_x}$)

With the saturated segments excluded, the remaining data satisfies the unsaturated dynamics (D). These data can now be used for least-squares estimation of $\theta_{NO_x}$.

### 3a. Regression form

The compact unsaturated dynamics $\bar{\eta}_F(k{+}1) = \bar{\eta}_F(k) + \phi_{NO_x}^T(k)\,\theta_{NO_x}$ (Section 4a of the derivation) rearrange into a standard linear regression by moving $\bar{\eta}_F(k)$ to the left:

$$
y_{NO_x}(k) = \phi_{NO_x}^T(k)\,\theta_{NO_x} \tag{R}
$$

where $y_{NO_x}(k) = \bar{\eta}_F(k{+}1) - \bar{\eta}_F(k)$. The regressor $\phi_{NO_x}$ depends only on measured quantities and is constructed to avoid unnecessary multiplication/division of the same signals (which would amplify noise).

### 3b. Estimation on unsaturated segments

$\theta_{NO_x}$ is estimated by least squares on the unsaturated segments from Step 2. Under Gaussian model structure error, least squares is the MLE. These are the same parameters that appear in the guard condition $G(k) = \phi_g^T(k)\,\theta_{NO_x}$ (Section 4c of the derivation), so once estimated, the guard is fully determined.

### 3c. Estimation results

With both $\hat{\theta}_{sat}$ and $\hat{\theta}_{NO_x}$, the full switched model (Section 4d) is simulated from initial conditions and inputs alone on the training data. Goodness of fit: $\%\text{fit} = 100 \times \bigl(1 - \frac{\|Y - \hat{Y}\|}{\|Y - \text{mean}(Y)\|}\bigr)$.

| Age | Test | Switched model | CSTR baseline |
|-----|------|---------------|---------------|
| Degreened | RMC | 71.9% | 20.0% |
| Aged | RMC | 69.9% | 27.7% |
| Degreened | hot-FTP | 61.9% | 7.1% |
| Aged | hot-FTP | 63.5% | 7.4% |
| Degreened | cold-FTP | 52.9% | 23.2% |
| Aged | cold-FTP | 60.3% | 23.2% |

The CSTR baseline column shows the linearized single-CSTR model fit on the same data for comparison.

Simulation overlays (training data):

| Degreened | Aged |
|-----------|------|
| ![DG cold-FTP](../Thesis/figs/4-NOx_mdl/3_figs/eta_sim_dg_cftp.png) | ![Aged cold-FTP](../Thesis/figs/4-NOx_mdl/3_figs/eta_sim_aged_cftp.png) |
| ![DG hot-FTP](../Thesis/figs/4-NOx_mdl/3_figs/eta_sim_dg_hftp.png) | ![Aged hot-FTP](../Thesis/figs/4-NOx_mdl/3_figs/eta_sim_aged_hftp.png) |
| ![DG RMC](../Thesis/figs/4-NOx_mdl/3_figs/eta_sim_dg_rmc.png) | ![Aged RMC](../Thesis/figs/4-NOx_mdl/3_figs/eta_sim_aged_rmc.png) |

---

## 4. Cross-validation

The parameter estimates from the "Degreened" training data are applied to three additional degreened catalyst runs (DG_1, DG_2, DG_3) to test generalization:

| Age | Test | %fit |
|-----|------|------|
| Degreened (train) | RMC | 71.86 |
| DG_1 | RMC | 67.12 |
| DG_2 | RMC | 54.71 |
| DG_3 | RMC | 67.01 |
| Degreened (train) | hot-FTP | 61.86 |
| DG_1 | hot-FTP | 54.87 |
| DG_2 | hot-FTP | 59.44 |
| DG_3 | hot-FTP | 50.17 |
| Degreened (train) | cold-FTP | 52.87 |
| DG_1 | cold-FTP | 39.41 |
| DG_2 | cold-FTP | 41.52 |
| DG_3 | cold-FTP | 42.17 |

The model generalizes across runs, with cross-validation %fit values within ~10-15 percentage points of the training fit.

---

*Compiled for the SCR-dynamics JDSMC publication.*
