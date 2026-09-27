# Conclusions and Application to Aging Diagnostics

> Continuation of [NOx_Dynamic_Model_Derivation.md](NOx_Dynamic_Model_Derivation.md) and [System_Identification.md](System_Identification.md).  
> Equations (S), (N), (D), (SW), $\theta_{sat}$, $\theta_{NO_x}$, $\phi_{sat}$, $\phi_2$, $\alpha_{sat}$, and the LP/LS formulations are defined there.  
> Source: Thesis Chapters 5 & 6.

---

## 1. From model parameters to aging diagnostics

### 1a. Physical basis: aging reduces adsorption site concentration

Catalyst aging (hydrothermal degradation) reduces the concentration of viable adsorption sites $\Gamma$ on the catalyst surface. Since $\Gamma$ appears directly in the saturated model ceiling:

$$
\eta_{sat}(k+1) = \frac{u_1(k)}{F(k)}\,\tau_0\,k_{scr}(T(k))\,\Gamma(T(k)) = \phi_{sat}^T(k)\,\theta_{sat}
$$

a decrease in $\Gamma$ lowers $\eta_{sat}$, i.e., the maximum achievable NOx reduction at any given operating condition drops. This effect is captured entirely by $\theta_{sat}$.

### 1b. The aging factor $\alpha_{sat}$

Normalizing $\eta_{sat}$ by inlet NOx and flow rate isolates a quantity that depends only on temperature and the product $k_{scr}\,\Gamma$:

$$
\alpha_{sat}(k) = \frac{F(k)}{u_1(k)}\,\eta_{sat}(k) = \tau_0\,k_{scr}(T(k))\,\Gamma(T(k)) = \phi_2^T(T(k))\,\theta_{sat}
$$

When plotted as a function of temperature, $\alpha_{sat}$ is consistently and statistically lower for an aged catalyst than for a degreened catalyst in both Hot-FTP and RMC tests, confirming that aging reduces $\Gamma$.

| Degreened (FTIR) | Aged (FTIR) |
|---|---|
| ![alpha_sat FTIR](../Thesis/figs/5-sat_detect/3_parm_ID/aging_factor_ssd.png) | ![alpha_sat NOx sensor](../Thesis/figs/5-sat_detect/3_parm_ID/aging_factor_iod.png) |

![Predicted eta_sat for aged vs. degreened](../Thesis/figs/5-sat_detect/3_parm_ID/AgingDiff.png)

### 1c. Why the saturated model is preferred for diagnostics

The saturated model parameters are preferred over the unsaturated model parameters for aging diagnostics for three reasons:

1. **No backward differencing.** The saturated model regression $\eta_{sat}(k) = \phi_{sat}^T(k)\,\theta_{sat}$ does not involve backward differences of the output, unlike the unsaturated regression (D) which computes $\bar{\eta}_F(k{+}1) - \bar{\eta}_F(k)$ and amplifies noise.
2. **Fewer parameters.** $\theta_{sat}$ has 3 parameters versus 9 for $\theta_{NO_x}$, requiring a less stringent persistence of excitation condition and reducing overfitting risk.
3. **Robustness to cross-sensitivity.** The LP estimation uses only data near saturation. At saturation, tailpipe NH$_3$ slip is minimal, so NOx sensor cross-sensitivity error is small compared to the unsaturated regime.

---

## 2. Hypothesis testing framework (Wald test)

### 2a. Formulation

Aging detection is posed as a composite hypothesis test:

$$
\mathcal{H}_0: \hat{\theta}_{sat} = \theta_{sat}^{dg} \quad \text{(degreened)}, \qquad
\mathcal{H}_1: \hat{\theta}_{sat} \neq \theta_{sat}^{dg} \quad \text{(aged)}
$$

where $\theta_{sat}^{dg}$ is the known parameter estimate for a degreened catalyst under the same input conditions (temperature, flow rate, urea dosing ranges).

The LP estimate $\hat{\theta}_{sat}$ is the MLE under a half-normal distribution for the model structure error $\varepsilon_\eta(k) = \eta_{sat}(k) - \eta(k) \ge 0$. By asymptotic MLE theory:

$$
\hat{\theta}_{sat} \sim \mathcal{N}\bigl(\theta_{sat},\; I^{-1}(\theta_{sat})\bigr), \qquad
I(\theta_{sat}) = \frac{1}{\sigma^2}\,\Phi_{sat}^T\,\Phi_{sat}
$$

The generalized likelihood ratio test reduces to the **Wald test statistic**:

$$
T_w = \bigl(\hat{\theta}_{sat} - \theta_{sat}^{dg}\bigr)^T\,I(\hat{\theta}_{sat})\,\bigl(\hat{\theta}_{sat} - \theta_{sat}^{dg}\bigr)
$$

Decide $\mathcal{H}_1$ (aged) if $T_w > \beta$.

### 2b. Distribution of the test statistic

Under $\mathcal{H}_0$, $T_w \sim \chi^2_k$ (central, $k = \dim(\theta_{sat}) = 3$). Under $\mathcal{H}_1$, $T_w \sim \chi^2_k(\lambda)$ (non-central), where:

$$
\lambda = \bigl(\theta_{sat}^{ag} - \theta_{sat}^{dg}\bigr)^T\,I(\theta_{sat}^{ag})\,\bigl(\theta_{sat}^{ag} - \theta_{sat}^{dg}\bigr)
$$

The non-centrality parameter $\lambda$ grows with the magnitude of aging-induced parameter shift. Performance measures:

$$
P_{FA} = 1 - F_{\chi^2_k}(\beta), \qquad
P_D = 1 - F_{\chi^2_k(\lambda)}(\beta), \qquad
P_{MD} = 1 - P_D
$$

The threshold $\beta$ can be chosen via the Neyman-Pearson criterion to achieve a desired false alarm probability and maximize detection probability.

### 2c. Unsaturated model detector (with nuisance parameters)

Among the unsaturated parameters $\theta_{NO_x} = [\theta_{\eta_{ads}}^T \; \theta_{\eta_{od}}^T \; \theta_{\eta_{scr}}^T \; \theta_\Gamma^T]^T$, only $\theta_\Gamma$ depends on $\Gamma$ and is thus sensitive to aging. The remaining parameters are nuisance parameters $\theta_\nu = [\theta_{\eta_{ads}}^T \; \theta_{\eta_{od}}^T \; \theta_{\eta_{scr}}^T]^T$ that must be accounted for:

$$
T_w = \bigl(\hat{\theta}_\Gamma - \theta_\Gamma^{dg}\bigr)^T\,\bigl[I^{-1}(\hat{\theta})\bigr]_{\theta_\Gamma\theta_\Gamma}^{-1}\,\bigl(\hat{\theta}_\Gamma - \theta_\Gamma^{dg}\bigr)
$$

where $\bigl[I^{-1}(\hat{\theta})\bigr]_{\theta_\Gamma\theta_\Gamma}^{-1} = I_{\theta_\Gamma\theta_\Gamma} - I_{\theta_\Gamma\theta_\nu}\,I_{\theta_\nu\theta_\nu}^{-1}\,I_{\theta_\nu\theta_\Gamma}$ is the Schur complement, accounting for the information loss due to the nuisance parameters. Under Gaussian model error, LS is the MLE for the unsaturated regression.

---

## 3. Diagnostic results

### 3a. Test-cell data: parameter estimates

Parameter estimates from test-cell experiments with FTIR measurements:

| Test | $\hat{\theta}_1$ | $\hat{\theta}_2$ | $\hat{\theta}_3$ | $\sigma$ |
|------|---|---|---|---|
| DG RMC 1 | -0.63 | 0.55 | 38.74 | 1.14 |
| DG RMC 2 | -0.82 | 0.69 | 38.61 | 1.08 |
| DG RMC 3 | -0.80 | 1.12 | 38.78 | 1.10 |
| Aged RMC | -1.31 | 1.89 | 36.94 | 1.01 |
| DG Hot FTP 1 | -2.14 | -5.80 | 39.53 | 1.35 |
| DG Hot FTP 2 | -3.45 | -7.99 | 38.05 | 1.30 |
| DG Hot FTP 3 | -3.37 | -6.42 | 38.88 | 1.33 |
| Aged Hot FTP | -4.54 | -8.01 | 37.00 | 1.31 |

Parameters show consistency across repeated tests on the same catalyst but vary between RMC and Hot-FTP due to differences in temperature ranges (RMC: 240-360°C; Hot-FTP: 200-300°C). The aged catalyst consistently shows a lower $\hat{\theta}_3$ (the constant term in the Chebyshev expansion), reflecting reduced $k_{scr}\,\Gamma$.

Parameter estimates using NOx sensor data:

| Test | $\hat{\theta}_1$ | $\hat{\theta}_2$ | $\hat{\theta}_3$ | $\sigma$ |
|------|---|---|---|---|
| DG RMC 1 | -0.70 | 0.56 | 38.48 | 1.10 |
| DG RMC 2 | -1.18 | 0.93 | 38.43 | 1.07 |
| DG RMC 3 | -0.54 | 0.62 | 38.99 | 1.06 |
| Aged RMC | -1.03 | 1.26 | 37.80 | 1.00 |
| DG Hot FTP 1 | -2.30 | -6.11 | 37.08 | 1.28 |
| DG Hot FTP 2 | -2.71 | -6.13 | 36.22 | 1.27 |
| DG Hot FTP 3 | -3.92 | -6.17 | 35.81 | 1.30 |
| Aged Hot FTP | -6.42 | -13.06 | 34.73 | 1.34 |

NOx sensor cross-sensitivity introduces bias into the parameter estimates (the portion of the cross-sensitivity error $\varepsilon_\chi(T)$ that projects onto the Chebyshev basis is absorbed into $\theta_{sat}$). This reduces the separation between aged and degreened estimates but preserves the overall aging trend.

### 3b. Test-cell data: Wald test results (saturated detector)

| Test | $T_w$ (FTIR) | $T_w$ (NOx sensor) |
|------|---|---|
| DG RMC 1 | 0.26 | 5.54 |
| DG RMC 2 | 0.00 | 0.00 |
| DG RMC 3 | 6.71 | 5.46 |
| Aged RMC | **47.16** | **17.68** |
| DG Hot FTP 1 | 0.49 | 0.81 |
| DG Hot FTP 2 | 3.26 | 1.98 |
| DG Hot FTP 3 | 0.00 | 0.00 |
| Aged Hot FTP | **4.84** | **11.44** |

Thresholds: $\beta_{RMC} = 10.0$, $\beta_{hFTP} = 4.0$. The aged catalyst exceeds the threshold in all cases. The non-central $\chi^2$ distribution under aging is different for RMC and FTP due to the differences in the ranges of inputs experienced during the tests.

### 3c. Test-cell data: Wald test results (unsaturated detector)

| Test | $T_w$ (FTIR) | $T_w$ (NOx sensor) |
|------|---|---|
| DG Hot FTP 1 | 2.95 | 2.72 |
| DG Hot FTP 2 | 13.13 | 8.81 |
| DG Hot FTP 3 | 0.00 | 0.00 |
| Aged Hot FTP | **164.14** | **165.33** |
| DG RMC 1 | 42.55 | 34.11 |
| DG RMC 2 | 0.00 | 0.00 |
| DG RMC 3 | 30.11 | 28.67 |
| Aged RMC | **119.61** | **111.56** |

The unsaturated detector works well on test-cell data (large separation for aged), but this is partly an artifact of the repeatable drive cycles producing uniform cross-sensitivity error profiles.

### 3d. Truck data: aging detection results

Road data from four long-haul trucks, collected at two times separated by years of operation. Drive segments are partitioned from engine-start to engine-stop; only segments with sufficient near-saturation samples ($N_{sat} \ge 150$, or $120$ for Truck C) are used. Truck threshold: $\beta_{truck} = 250$.

**Summary:**

| Truck | Total Segments | Detection | False Alarm | Missed Detection |
|-------|---|---|---|---|
| Truck A | 8 | 4 | 0 | 4 |
| Truck B | 11 | 9 | 1 | 1 |
| Truck C | 9 | 7 | 1 | 1 |
| Truck D | 8 | 6 | 2 | 0 |

Trucks B, C, D show high detection rates with low false-alarm and missed-detection rates. Truck A has lower performance, attributed to relatively lower mileage difference between the two data collection periods.

**Detailed segment results:**

**Truck A:**

| Segment | $T_w$ | $N_{sat}$ | $N_{sat}/N$ (%) |
|---------|---|---|---|
| A_15_0 | 0.00 | 797 | 6.70 |
| A_15_1 | 28.95 | 305 | 2.75 |
| A_15_2 | 58.50 | 276 | 5.93 |
| A_15_3 | 89.56 | 271 | 12.77 |
| A_17_0 | 65.59 | 414 | 2.97 |
| A_17_1 | 127.79 | 614 | 6.76 |
| A_17_2 | 12.13 | 330 | 8.44 |
| A_17_3 | 61.76 | 336 | 9.73 |

**Truck B:**

| Segment | $T_w$ | $N_{sat}$ | $N_{sat}/N$ (%) |
|---------|---|---|---|
| B_15_0 | 672.43 | 167 | 1.26 |
| B_15_1 | 106.66 | 189 | 1.81 |
| B_15_3 | 237.79 | 159 | 2.95 |
| B_15_4 | 0.00 | 226 | 6.03 |
| B_15_6 | 84.73 | 179 | 5.75 |
| B_18_0 | 543.55 | 231 | 2.10 |
| B_18_1 | 1050.24 | 308 | 3.00 |
| B_18_4 | 22.48 | 195 | 6.14 |
| B_18_6 | 7150.19 | 378 | 12.73 |
| B_18_7 | 760.89 | 224 | 8.14 |
| B_18_8 | 342.67 | 246 | 9.04 |

**Truck C:**

| Segment | $T_w$ | $N_{sat}$ | $N_{sat}/N$ (%) |
|---------|---|---|---|
| C_15_1 | 0.00 | 186 | 1.30 |
| C_15_4 | 113.86 | 143 | 3.84 |
| C_15_5 | 4340.08 | 124 | 4.16 |
| C_17_0 | 700.88 | 186 | 1.24 |
| C_17_2 | 3993.85 | 446 | 7.23 |
| C_17_3 | 1693.49 | 199 | 4.00 |
| C_17_4 | 186.64 | 120 | 2.81 |
| C_17_7 | 1878.52 | 129 | 4.73 |
| C_17_9 | 7500.10 | 218 | 11.63 |

**Truck D:**

| Segment | $T_w$ | $N_{sat}$ | $N_{sat}/N$ (%) |
|---------|---|---|---|
| D_15_0 | 0.00 | 214 | 1.29 |
| D_15_2 | 98.99 | 160 | 2.95 |
| D_15_4 | 1546.55 | 213 | 4.90 |
| D_15_6 | 935.05 | 152 | 8.86 |
| D_16_1 | 423.20 | 246 | 3.67 |
| D_16_2 | 1022.06 | 199 | 5.22 |
| D_16_3 | 697.78 | 170 | 5.36 |
| D_16_5 | 313.62 | 200 | 9.31 |

There is a strong correlation between detection of aging and the number of samples close to catalyst saturation in the drive segment. In truck data, fewer than 10% of samples are typically near saturation, unlike test-cell data where the fraction is much higher.

The unsaturated detector does not generalize reliably to truck data due to: (i) the 9-parameter model overfitting noise without sufficient excitation, (ii) backward differencing amplifying sensor noise, and (iii) varying cross-sensitivity error profiles across non-repeatable real-world drive cycles.

---

## 4. Conclusions

### 4a. Key contributions

1. **Discrete switched nonlinear model for SCR-ASC dynamics.** Derived from molar conservation (not CSTR linearization), the model incorporates the interplay between residence time and sensor sampling time, uses Chebyshev polynomial bases for linear-in-parameter structure, and models saturation/desaturation via a switched structure with a causal switching condition.

2. **Identifiable parametrization and convex estimation.** The unmeasured surface ammonia state $\sigma$ is eliminated algebraically via $\eta(k)$, yielding a recursive model linear in parameters. The saturated parameters are estimated via LP; unsaturated parameters via LS. Both are convex problems.

3. **Model validation.** The switched model achieves goodness-of-fit exceeding 50% on all test-cell cycles (RMC, Hot-FTP, Cold-FTP) for both degreened and aged catalysts, consistently outperforming the linearized single-CSTR model (typically below 30% on the same data).

4. **Statistical aging detection framework.** A Wald-test detector, derived within an MLE framework using the asymptotic properties of the Fisher information matrix, provides:
   - analytical distributions for the test statistic under both hypotheses ($\chi^2$ / non-central $\chi^2$);
   - threshold selection via the Neyman-Pearson criterion;
   - performance measures ($P_{FA}$, $P_D$, $P_{MD}$).

5. **Experimental validation.** The saturated-model detector reliably distinguishes aged from degreened catalysts on both test-cell data and real-world truck data (four long-haul trucks), despite NOx sensor cross-sensitivity and variable operating conditions.

### 4b. Limitations

1. **Reduced reaction set.** Only three reactions are retained (standard SCR, lumped AMOX, NH$_3$ adsorption/desorption). Fast/slow SCR, NO oxidation, and side reactions are omitted. The lumped AMOX with 100% N$_2$ selectivity is valid only above ~225°C.

2. **Temperature-range dependence.** Chebyshev polynomial approximations tie parameter estimates to the specific temperature window of the estimation data. A single parameter set does not extrapolate across widely different thermal profiles (e.g., RMC vs. Hot-FTP).

3. **Dependence on catalyst saturation data.** The detector requires sufficient near-saturation samples ($N_{sat} \ge N_{sat}^{min}$). In real-world driving, fewer than 10% of samples are typically near saturation, limiting the duty cycles over which the detector is operational.

4. **NOx sensor cross-sensitivity.** Cross-sensitivity reduces the separation between aged and degreened parameter estimates, weakening discriminative power. It is treated as an unknown bounded disturbance rather than being explicitly corrected (tailpipe NH$_3$ measurements are unavailable on production trucks).

5. **Asymptotic assumptions.** The Wald test relies on large-sample MLE properties. For short drive segments or limited near-saturation data, the $\chi^2$ approximation may not hold, and $P_{FA}$/$P_D$ deviate from analytical values.

6. **Absence of ground truth.** Truck data validation is limited by the absence of ground-truth aging levels; only mileage difference serves as a proxy. The relationship between $\lambda$ and true aging severity cannot be precisely calibrated.

7. **Threshold calibration.** Truck data thresholds were selected empirically from a limited dataset rather than via the Neyman-Pearson criterion, which requires sufficient data to characterize both hypothesis distributions.

### 4c. Other applications and future work

1. **Urea dosing control.** The switched nonlinear model structure can inform more efficient urea dosing control policies, e.g., by maximizing the duration of safely-maintained catalyst saturation while respecting ammonia slip and hardware constraints.

2. **SETAR framework.** The Self-Excited Threshold Autoregressive (SETAR) structure offers a promising avenue for capturing regime-switching behavior in SCR-ASC dynamics more generally.

3. **Wider temperature validity.** Alternative formulations for the temperature dependence of reaction and transport parameters could achieve a wider temperature validity range and unique identifiability across operating conditions.

4. **Extended reaction pathways.** Incorporating fast/slow SCR reactions and explicit ASC oxidation dynamics into the switched framework, especially where inter-brick or tailpipe ammonia sensors are available, would improve model fidelity at lower temperatures.

5. **Cross-sensitivity correction.** Explicit estimation and correction of the NOx sensor cross-sensitivity to ammonia would sharpen the separation between aged and degreened parameter estimates.

6. **Fleet-scale validation.** A broader validation campaign using a larger fleet with controlled, documented aging levels would enable rigorous Neyman-Pearson threshold calibration and quantify how $\lambda$ scales with true degradation severity.

7. **ECU integration.** Integration of the detector into an engine control unit prototype, with real-time parameter estimation running alongside the existing urea dosing controller, would assess computational feasibility, memory requirements, and latency.

---

*Compiled for the SCR-dynamics JDSMC publication.*
