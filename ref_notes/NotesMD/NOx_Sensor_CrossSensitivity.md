# NOx Sensor Cross-Sensitivity

> Consequences of NOx sensor cross-sensitivity to ammonia for model estimation and diagnostics.  
> Builds on [NOx_Dynamic_Model_Derivation.md](NOx_Dynamic_Model_Derivation.md) and [System_Identification.md](System_Identification.md).  
> Source: Thesis Chapters 2 (3-cross.tex), 5 (2-param_dist.tex, 3-testcell_aging.tex, 5-unconstrained_results.tex, 6-detector_performance.tex).

---

## 1. The cross-sensitivity error

Commercial NOx sensors are cross-sensitive to tailpipe ammonia. The measurement error is:

$$
\varepsilon_\chi(k) = \chi(T(k)) \cdot [\text{NH}_3]^{out}(k) \ge 0
$$

where $\chi(T)$ is a temperature-dependent cross-sensitivity factor. The error is **non-negative and directional**: the sensor always over-reports tailpipe NOx when ammonia is present.

The measured quantities are:

$$
y_1(k) = x_1(k) + \varepsilon_\chi(k), \qquad \eta_y(k) = u_1(k{-}1) - y_1(k) = \eta(k) - \varepsilon_\chi(k) \le \eta(k)
$$

The sensor-derived NOx reduction $\eta_y$ **under-estimates** the true reduction $\eta$.

![NOx sensor cross-sensitivity to tailpipe ammonia](../Thesis/figs/2-data/cross_aged_hftp.png)

---

## 2. Consequences for the saturated model ($\theta_{sat}$)

### 2a. Bounding condition is preserved

The bounding condition $\eta_{sat}(k) \ge \eta(k)$ (Section 1a of System_Identification.md) remains valid when using sensor data:

$$
\eta_{sat}(k) \ge \eta(k) = \eta_y(k) + \varepsilon_\chi(k) \ge \eta_y(k)
$$

Since $\varepsilon_\chi \ge 0$, the measured reduction $\eta_y$ is always below the true reduction, so the LP for $\theta_{sat}$ is still feasible and produces a valid upper bound.

### 2b. Parameter bias

The cross-sensitivity error has a component with the same temperature dependence as the model:

$$
\varepsilon_\chi(T) = \phi_{sat}^T(T)\,\theta_{\varepsilon_\chi} + \Delta\varepsilon_\chi
$$

The component $\theta_{\varepsilon_\chi}$ becomes an additive bias in the estimated $\theta_{sat}$. This reduces the separation between aged and degreened parameter estimates but preserves the directional trend (aged parameters are consistently lower).

### 2c. Saturated detector robustness

The saturated-model-based aging detector is **more robust** to cross-sensitivity than the unsaturated detector because it uses only data points near saturation. Near saturation, urea dosing is high and the catalyst surface is fully occupied, so the cross-sensitivity error has a more predictable structure.

---

## 3. Consequences for the unsaturated model ($\theta_{NO_x}$)

The unsaturated regression $y_{NO_x} = \phi_{NO_x}^T\,\theta_{NO_x}$ (R) involves a backward difference in $\bar{\eta}_F$, which amplifies measurement noise. Cross-sensitivity compounds this problem.

**Test-cell data:** The Wald test for aging detection using $\theta_{NO_x}$ works on test-cell data even with NOx sensor cross-sensitivity. This is because the drive cycles are repeatable, so $\varepsilon_\chi(k)$ has a similar time response across runs. The cross-sensitivity effect partially cancels when comparing aged vs. degreened parameter estimates from the same test protocol.

**Truck data:** The test is inconclusive. Three limitations are identified:
1. Cross-sensitivity time responses vary across drive segments (non-repeatable real-world driving).
2. The unsaturated model has 9 parameters (risk of overfitting with limited data).
3. Backward difference in the regressor amplifies both sensor noise and cross-sensitivity.

---

## 4. Cross-sensitivity estimation

### 4a. FTIR bias correction

The FTIR reference sensor has a time-varying bias $b(t) = b_1\,t + b_0$ (drift + offset). Estimated from the starting and tail segments of the data. Result: drift $b_1$ is negligible compared to bias $b_0$.

![Sensor bias correction for RMC cycles](../Thesis/figs/2-data/chi_est/dg_rmc_NOx_bias.eps)

### 4b. Constant $\chi$ (temperature-independent)

For the RMC temperature range, $\chi$ is approximately constant. Regression:

$$
\underbrace{y_1 - x_1}_{\mathbf{y}} = \underbrace{\begin{bmatrix} x_2 & -1 \end{bmatrix}}_{\phi^T} \underbrace{\begin{bmatrix} \chi \\ \chi\,b_{th} + b_0 \end{bmatrix}}_{\theta}
$$

where $b_{th}$ is the ammonia threshold below which cross-sensitivity is negligible.

### 4c. Temperature-dependent $\chi(T) = aT - b_T$

Regression with temperature interaction:

$$
y_1 - x_1 = \begin{bmatrix} T\,x_2 & -T & -x_2 & 1 \end{bmatrix} \begin{bmatrix} a \\ a\,b_{th} \\ b_T \\ b_T\,b_{th} - b_0 \end{bmatrix}
$$

Temperature modeling reduces prediction error, but the improvement over the constant model is not significant.

### 4d. Decision

Tailpipe ammonia measurements ($x_2$) are unavailable on commercial vehicles. Therefore, **cross-sensitivity is treated as an unknown bounded disturbance** rather than being estimated and corrected.

---

## 5. Summary of cross-sensitivity effects

| Aspect | Effect | Mitigation |
|--------|--------|-----------|
| Measured $\eta_y$ | Under-estimates true $\eta$ | Non-negative error preserves bounding condition |
| $\theta_{sat}$ estimation (LP) | Additive bias from $\varepsilon_\chi$ temperature component | Aging trend preserved; saturated detector robust |
| $\theta_{NO_x}$ estimation (LS) | Backward difference amplifies noise + cross-sensitivity | Works on test-cell (repeatable cycles); fails on truck data |
| Aging detection | Reduced aged/degreened separation | Saturated model detector preferred |
| Estimation of $\chi$ itself | Feasible with FTIR data; $\chi(T)$ marginally better than constant | Not feasible on trucks (no $\text{NH}_3$ sensor) |

---

*Compiled for the SCR-dynamics JDSMC publication.*
