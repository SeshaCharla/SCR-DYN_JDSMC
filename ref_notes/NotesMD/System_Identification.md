# System Identification for the NOx Dynamic Model

> From the SETNARX model to parameter estimates.  
> Source: Thesis Chapters 2, 4, 5.

---

The switched model has two parameter vectors: $\theta_{sat}$ (saturation ceiling) and $\theta_{NO_x}$ (unsaturated dynamics and guard condition). The goal is to estimate both from available measurement data ($u_1, x_1, u_2, F, T$).

We address this by first estimating $\theta_{sat}$ using a linear programming approach that does not require mode knowledge, then using the estimated saturation ceiling to identify data segments near saturation, and finally estimating $\theta_{NO_x}$ from the remaining unsaturated data.

---

## 1. Saturated model parameters ($\theta_{sat}$)

### 1a. Regression form

Under catalyst saturation ($\sigma = \Gamma$), the NOx reduction depends only on the current operating point (not on the previous state or urea dosing):

$$
\eta_{sat}(k{+}1) = \frac{u_1(k)}{F(k)}\,\tau_0\,k_{scr}(T(k))\,\Gamma(T(k))
$$

Parametrizing $k_{scr}\,\Gamma$ with the quadratic Chebyshev basis $\phi_2$:

$$
\eta_{sat}(k{+}1) = \phi_{sat}^T(k)\,\theta_{sat}, \qquad \phi_{sat}(k) = \frac{u_1(k)}{F(k)}\,\phi_2(T(k))
$$

where $\theta_{sat} = \Gamma\,\tau_0\,k_{s2v}\,\theta_{scr}$. The flow-normalized quantity $\alpha_{sat}(k) = \frac{F(k)}{u_1(k)}\,\eta_{sat}(k) = \phi_2^T(T)\,\theta_{sat}$ is purely a function of temperature.

> **Thesis plots:** $\alpha_{sat}$ vs temperature for aged/degreened catalysts shows a consistent and statistically significant decrease with aging (Fig. alpha_sat.png, alpha_sat_1.png). This supports the hypothesis that aging reduces $\Gamma$.

### 1b. Bounding condition

For any time step $k$, the saturated ceiling bounds the actual NOx reduction:

$$
\eta_{sat}(k) \ge \eta(k) \qquad \forall\, k
$$

This holds because $\sigma = \Gamma$ always produces at least as much reduction as any $\sigma \le \Gamma$. The bound is mode-independent, so all data can be used.

### 1c. LP formulation

The tightest upper bound minimizes the total area under the saturated curve:

$$
\hat{\theta}_{sat} = \arg\min_{\theta_{sat}} \;\mathbf{1}^T\,\Phi_{sat}\,\theta_{sat} \qquad \text{s.t.} \quad \Phi_{sat}\,\theta_{sat} \succeq H
$$

where $\Phi_{sat} = [\phi_{sat}(1), \ldots, \phi_{sat}(N{-}1)]^T$ and $H = [\eta(2), \ldots, \eta(N)]^T$.

> **Thesis plots:** Bounded $\eta$ plots for degreened/aged x RMC/hot-FTP/cold-FTP (eta_bounds_*.png). The saturated ceiling forms a tight envelope over the actual $\eta$ data. The maximum NOx reduction is inversely related to flow rate (higher residence time allows more reduction).

### 1d. MLE interpretation (half-normal error)

The model structure error $\varepsilon_\eta(k) = \eta_{sat}(k) - \eta(k) \ge 0$ is non-negative by construction. Modelling $\varepsilon_\eta$ as half-normal with scale $\sigma$:

$$
p(\varepsilon_\eta;\,\sigma) = \frac{\sqrt{2}}{\sigma\sqrt{\pi}}\,\exp\!\Bigl(-\frac{\varepsilon_\eta^2}{2\sigma^2}\Bigr), \qquad \varepsilon_\eta \ge 0
$$

> **Thesis plot:** Distribution of $\varepsilon_\eta$ (eta_dist.png) validates the half-normal assumption: non-negative with a long tail.

Maximizing the log-likelihood:

$$
L(\theta_{sat}) = \frac{N}{2}\ln\frac{2}{\sigma^2\pi} - \frac{1}{2\sigma^2}\sum_{k=1}^{N}\varepsilon_\eta^2(k), \qquad \varepsilon_\eta(k) \ge 0
$$

is equivalent to minimizing $\sum \varepsilon_\eta^2$ s.t. $\varepsilon_\eta \ge 0$ (constrained QP). By norm equivalence ($l_2 \to l_1$) and the non-negativity constraint, the QP relaxes to the LP.

**Result:** The LP solution is the MLE of $\theta_{sat}$ under half-normal model structure error.

### 1e. Parameter distribution

By asymptotic MLE theory (large $N$):

$$
\hat{\theta}_{sat} \sim \mathcal{N}\bigl(\theta_{sat},\; I^{-1}(\theta_{sat})\bigr), \qquad I(\theta_{sat}) = \frac{1}{\sigma^2}\,\Phi_{sat}^T\,\Phi_{sat}
$$

Prediction variance: $\text{Var}[\eta_{sat}(k)] = \phi_{sat}^T(k)\,I^{-1}(\theta_{sat})\,\phi_{sat}(k)$.

> **Thesis tables:** $\theta_{sat}$ estimates for degreened/aged x RMC/hot-FTP/cold-FTP with two temperature zones (Tab. sat_parm_est). Parameter values are directionally different for aged vs. degreened.

### 1f. NOx sensor cross-sensitivity

Commercial NOx sensors have cross-sensitivity to tailpipe ammonia: $y_1(k) = x_1(k) + \chi(T)\,[\text{NH}_3]^{out}(k)$, introducing a non-negative error $\varepsilon_\chi \ge 0$. The measured NOx reduction underestimates the true value: $\eta_y(k) = \eta(k) - \varepsilon_\chi(k) \le \eta(k)$. Since $\eta_{sat} \ge \eta \ge \eta_y$, the bounding condition and LP remain valid when using sensor data.

---

## 2. Segment identification using $\theta_{sat}$

With $\hat{\theta}_{sat}$ estimated, the saturated model response $\hat{\eta}_{sat}(k) = \phi_{sat}^T(k)\,\hat{\theta}_{sat}$ is computed for each time step. A data point is classified as near-saturated if the prediction error is small:

$$
|\hat{\eta}_{sat}(k) - \eta(k)| \le \varepsilon_{sat} \implies \text{saturated operation}
$$

The threshold $\varepsilon_{sat}$ is chosen from the variance of the prediction error distribution ($\varepsilon_{sat} = 2.5 \times 10^{-3}\;\text{mol/m}^3$ for the test-cell dataset).

The saturated model response is independent of the previous state and urea dosing, depending only on $u_1$, $F$, and $T$. Saturation tends to occur when mass flow and urea dosing are high, consistent with physical intuition.

> **Thesis plot:** Catalyst mode detection for RMC data (SatDetect_dg_rmc_1_FTIR.png) showing identified saturated segments overlaid on the data.  
> **Thesis tables:** $N_{sat}/N$ fractions per test (Tab. sat_data_perc_ssd, sat_data_perc_iod). RMC tests have ~20% saturated data; hot-FTP has ~30%.

The unsaturated segments (the complement) are used for $\theta_{NO_x}$ estimation in the next step.

---

## 3. Unsaturated model parameters ($\theta_{NO_x}$)

### 3a. Regression form

The unsaturated dynamics (D) rearrange into a linear regression:

$$
y_{NO_x}(k) = \phi_{NO_x}^T(k)\,\theta_{NO_x} \tag{R}
$$

where:

$$
y_{NO_x}(k) = \frac{F(k)}{u_1(k)}\,\eta(k{+}1) - \frac{F(k{-}1)}{u_1(k{-}1)}\,\eta(k)
$$

$$
\phi_{NO_x}(k) = \begin{bmatrix}
\bigl(\frac{x_1(k)}{u_1(k{-}1)} - 1\bigr)\,u_2(k{-}1)\,\phi_1^T(k{-}1) \\[4pt]
\bigl(\frac{x_1(k)}{u_1(k{-}1)} - 1\bigr)\,F(k{-}1)\,\phi_1^T(k{-}1) \\[4pt]
-\eta(k)\,F(k{-}1)\,\phi_1^T(k) \\[4pt]
\frac{u_2(k{-}1)}{F(k{-}1)}\,\phi_2(k{-}1)
\end{bmatrix}, \qquad
\theta_{NO_x} = \begin{bmatrix} \theta_{\eta_{ads}} \\ \theta_{\eta_{od}} \\ \theta_{\eta_{scr}} \\ \theta_\Gamma \end{bmatrix}
$$

The regressor $\phi_{NO_x}$ is constructed to avoid unnecessary multiplication and division of the same signals, which would amplify noise.

### 3b. Estimation on unsaturated segments

$\theta_{NO_x}$ is estimated by least squares on the unsaturated segments identified in Step 2. Under Gaussian model structure error, least squares is the MLE. These are the same parameters that appear in the guard condition $G(k) = \phi_g^T(k)\,\theta_{NO_x}$, so once estimated, the guard condition is fully determined.

---

## 4. Validation

With both $\hat{\theta}_{sat}$ and $\hat{\theta}_{NO_x}$ estimated, the full switched model is simulated from initial conditions and inputs ($u_1, u_2, F, T$) alone. The guard condition selects the active mode at each step. Tailpipe NOx: $x_1(k{+}1) = u_1(k) - \eta(k{+}1)$.

Goodness of fit: $\%\text{fit} = 100 \times \bigl(1 - \frac{\|Y - \hat{Y}\|}{\|Y - \text{mean}(Y)\|}\bigr)$.

> **Thesis results (Tab. results):**
> | Age | Test | Switched model | CSTR model |
> |-----|------|---------------|------------|
> | Degreened | RMC | 71.9% | 20.0% |
> | Aged | RMC | 69.9% | 27.7% |
> | Degreened | hot-FTP | 61.9% | 7.1% |
> | Aged | hot-FTP | 63.5% | 7.4% |
> | Degreened | cold-FTP | 52.9% | 23.2% |
> | Aged | cold-FTP | 60.3% | 23.2% |
>
> The switched nonlinear model consistently fits >50% across all tests, significantly outperforming the linearized CSTR model.

> **Thesis plots:** Simulation overlays for all six test cases (eta_sim_*.png). Cross-validation on additional degreened data (cross_valid.png) confirms parameter generalization.

---

*Compiled for the SCR-dynamics JDSMC publication.*
