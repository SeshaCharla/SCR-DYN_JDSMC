# Practical Implementation Considerations

> Practical constraints, sensor effects, and preprocessing requirements for implementing the NOx dynamic model.  
> Builds on [NOx_Dynamic_Model_Derivation.md](NOx_Dynamic_Model_Derivation.md) and [System_Identification.md](System_Identification.md).  
> Source: Thesis Chapters 2, 3, 4, 5.

---

## 1. NOx sensor cross-sensitivity

Commercial NOx sensors are cross-sensitive to tailpipe ammonia:

$$
y_1(k) = x_1(k) + \varepsilon_\chi(k), \qquad \varepsilon_\chi(k) = \chi(T(k)) \cdot [\text{NH}_3]^{out}(k) \ge 0
$$

The error is non-negative and directional: the sensor always over-reports NOx when ammonia is present. The measured NOx reduction under-estimates the true value: $\eta_y(k) = \eta(k) - \varepsilon_\chi(k) \le \eta(k)$.

**Effect on $\theta_{sat}$ (LP):** The bounding condition $\eta_{sat} \ge \eta \ge \eta_y$ is preserved. LP remains feasible. Cross-sensitivity adds a bias to $\theta_{sat}$ (the component of $\varepsilon_\chi$ with the same temperature basis gets absorbed), reducing aged/degreened separation but preserving the directional trend.

**Effect on $\theta_{NO_x}$ (LS):** The backward difference in the unsaturated regressor amplifies both sensor noise and cross-sensitivity. Works on test-cell data (repeatable drive cycles give similar $\varepsilon_\chi$ profiles across runs) but is inconclusive on truck data.

**Decision:** Tailpipe NH$_3$ measurements are unavailable on commercial vehicles. Cross-sensitivity is therefore treated as an unknown bounded disturbance rather than being estimated and corrected.

---

## 2. Sensor bias and drift

FTIR reference sensors have a time-varying bias $b(t) = b_1\,t + b_0$. Estimated from the start and tail segments of the data. Finding: drift $b_1$ is negligible compared to offset $b_0$.

Truck data has **no FTIR and no NH$_3$ sensor**, so bias correction and cross-sensitivity estimation are infeasible on trucks.

---

## 3. Temperature range bounding

The polynomial rate-constant models have bounded validity. Data must be clipped to the operating range before estimation:

| Test | Temperature range |
|------|-------------------|
| RMC | 250--350 $°$C |
| hot-FTP, truck | 200--300 $°$C |

The linear Arrhenius approximation is valid $\pm 50°$C around $T_0 = 250°$C; the quadratic extends validity. Parameters do not generalize across thermal profiles. No extrapolation beyond the fitted window.

The temperature range is partitioned into two zones for $\theta_{sat}$ because $k_{scr} \cdot \Gamma$ has opposing monotonicity (k_{scr} increases, $\Gamma$ decreases with $T$), creating an inflection that a single polynomial cannot capture.

---

## 4. Flow rate constraints

The urea model $[\text{NH}_3]^{in} = \nu_u\,u_2/F$ has a reciprocal singularity at $F = 0$. Samples near zero flow must be excluded.

High flow creates a **virtual ceiling** on $\eta$ (reduced residence time limits NOx reduction) that is indistinguishable from catalyst saturation. This is a potential false positive for mode detection.

---

## 5. Sampling and residence time

| Platform | Sampling rate | Mode residence time |
|----------|---------------|---------------------|
| Test cell | 5 Hz ($t_s = 0.2$ s) | ~0.07 s |
| Truck | 1 Hz ($t_s = 1$ s) | ~0.07 s |

The model assumes an integer number of residence times per sample ($n = t_s/\tau$). Non-integer residuals are absorbed into model structure error. Assumption A9 ($T(k) \approx T(k{-}1)$) requires temperature to change slowly relative to the sampling period.

Test-cell data is decimated from 5 Hz to 1 Hz to match truck sampling. Upsampling truck data is not possible (unrecoverable high-frequency content).

---

## 6. Data preprocessing pipeline

1. **Unit conversion:** concentrations to $10^{-3}$ mol/m$^3$ (ppm $\to$ mol/m$^3$ requires gas temperature), $T$ offset by +200°C, $F$ scaled by $\times 10$ g/s.
2. **NaN filtering:** remove missing samples (common in truck data).
3. **Temperature clipping:** restrict to the operating window (Section 3).
4. **Decimation:** test-cell 5 Hz $\to$ 1 Hz.
5. **Smoothing:** non-causal 7th-order Chebyshev filter, 0--0.1 Hz passband. Truck gaps $< 1$ min bridged by linear interpolation before smoothing.
6. **Drive segment generation:** truck data split into contiguous intervals without breaks (not aligned to engine on/off cycles).
7. **Derived signals:** $\eta(k) = u_1(k{-}1) - x_1(k)$, $u_{2F} = u_2/F$, $u_{1F} = u_1/F$, $\bar{\eta}_F = \eta/u_{1F}(k{-}1)$. Note the one-sample delay convention in $\eta$.

---

## 7. Minimum data requirements

The LP for $\theta_{sat}$ **fails if the data lacks saturated regions**: the LP then fits the maximum observed $\sigma$ as if it were $\Gamma$, producing biased estimates.

For asymptotic MLE/Wald test validity, a minimum number of near-saturation samples is required:
- $N_{sat} \ge 150$ (120 for Truck C) per segment
- Segments with fewer saturated samples are pruned

More saturation data $\to$ better parameter estimates. RMC tests have ~20% saturated data; hot-FTP has ~30%.

---

## 8. Numerical issues

- **Backward-difference noise amplification:** the unsaturated regressor $\phi_{NO_x}$ contains a backward difference of the output, amplifying sensor noise. The saturated model avoids this entirely.
- **Signal ratio noise:** regressors are constructed to cancel unnecessary $u_1/F$ ratios (avoid multiplying and dividing the same measured signal).
- **Chebyshev conditioning:** mapping $T \to [-1,1]$ via Chebyshev polynomials instead of monomials prevents ill-conditioning in the temperature basis.
- **1/F singularity:** exclude $F \approx 0$ samples.
- **Max-min clipping:** simulation must enforce $\eta(k{+}1) = \max\{0, \min\{f_\sigma, f_\Gamma\}\}$ and also $x_1(k{+}1) = \max\{0, u_1(k) - \eta(k{+}1)\}$ to handle the case where potential reduction exceeds inlet NOx.

---

## 9. Chebyshev polynomial order

- $\phi_1$ (linear): individual rate constants ($k_{ads}$, $k_{od}$, $k_{scr}$), which are nearly linear in $T$ over the operating range.
- $\phi_2$ (quadratic): products involving $\Gamma$ ($k_{scr/ads} \cdot \Gamma$, $k_{scr} \cdot \Gamma$), which exhibit an inflection due to opposing monotonicity.

---

## 10. Drive cycle considerations

| Cycle | Characteristics | Saturated data |
|-------|----------------|----------------|
| RMC | Highest temperatures, repeatable | ~20% |
| hot-FTP | Moderate temperatures, repeatable | ~30% |
| cold-FTP | Low temperatures, repeatable | minimal |
| Truck | 200--300°C, non-repeatable segments | $\le 10\%$ |

Test-cell cycle repeatability means cross-sensitivity error profiles are similar across runs, which inflates unsaturated-model detector performance. This is an artifact not reproducible on truck data.

Truck segments are non-repeatable, unbounded length, and have no engine on/off alignment. No ground-truth aging measurement exists for truck data.

---

## 11. Practical thresholds

| Parameter | Value | Basis |
|-----------|-------|-------|
| $\varepsilon_{sat}$ | $2.5 \times 10^{-3}$ mol/m$^3$ | Prediction error variance |
| $\beta_{hFTP}$ | 4.00 | Empirical (Wald test threshold) |
| $\beta_{RMC}$ | 10.00 | Empirical |
| $\beta_{truck}$ | 250 | Empirical |
| $N_{sat}^{min}$ | 150 | Asymptotic MLE validity |

All thresholds are dataset-specific. Recalibrate empirically (bootstrap/Monte Carlo) for new fleets; asymptotic $\chi^2$ guarantees require large $N_{sat}$.

---

*Compiled for the SCR-dynamics JDSMC publication.*
