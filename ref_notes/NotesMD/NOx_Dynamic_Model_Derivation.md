# NOx Dynamic Model — Compact Derivation

> From reactions to the switched nonlinear NOx model (SETNARX).  
> Source: Thesis Chapters 3 & 4.

---

## Notation

| Symbol | Meaning |
|--------|---------|
| $[\cdot]$, $[\cdot]^{ads}$ | Volumetric (mol/m³) and surface (mol/m²) concentrations |
| $\{\cdot\}$ | Moles |
| $\Gamma$ | Total surface site concentration (mol/m²) |
| $k_{s2v} = A_{scr}/V$ | Surface-to-volume factor (m⁻¹) |
| $\tau = \tau_0/F$ | Residence time; $\tau_0 = V\rho_0$ |
| $t_s = n\tau$ | Sampling period ($n$ residence times per sample) |
| $u_1 = [\text{NO}_x]^{in}$, $x_1 = [\text{NO}_x]^{out}$ | Inlet / outlet NOx |
| $u_2$ | Urea injection rate (ml/s); $F$ = exhaust mass flow; $T$ = temperature |
| $k_{od} = k_{oxi} + k_{des}$ | Lumped oxidation + desorption rate constant |

---

## 1. Reactions and rate equations

The SCR-ASC system involves urea decomposition, ammonia adsorption/desorption, selective catalytic reduction of NOx, and ammonia oxidation (AMOX) on both SCR and ASC surfaces. For a tractable model, we retain only the three dominant reactions. The fast, slow, and standard SCR reactions are lumped into a single SCR reaction with an effective rate constant. The AMOX reactions on both SCR and ASC are lumped into a single oxidation reaction with 100% N$_2$ selectivity (valid for $T > 225\,°\text{C}$). The Eley–Rideal mechanism is assumed: gaseous NOx reacts with surface-adsorbed NH$_3$.

$$
4\,\text{NH}_3^{ads} + 4\,\text{NO} + \text{O}_2 \xrightarrow{k_{scr}} 4\,\text{N}_2 + 6\,\text{H}_2\text{O}, \qquad
4\,\text{NH}_3^{ads} + 3\,\text{O}_2 \xrightarrow{k_{oxi}} 2\,\text{N}_2 + 6\,\text{H}_2\text{O}, \qquad
\text{NH}_3 + \Theta_{free} \xrightleftharpoons[k_{des}]{k_{ads}} \text{NH}_3^{ads}
$$

The surface NH$_3$ rate follows directly from these three reactions — adsorption fills sites proportional to free-site concentration $(\Gamma - [\text{NH}_3]^{ads})$, while desorption, SCR, and oxidation each consume adsorbed NH$_3$:

$$
\dot{[\text{NH}_3]}^{ads} = \underbrace{k_{ads}\,[\text{NH}_3]^{in}\,(\Gamma - [\text{NH}_3]^{ads})}_{\text{adsorption}} - \underbrace{(k_{des} + k_{scr}\,[\text{NO}_x]^{in} + k_{oxi})\,[\text{NH}_3]^{ads}}_{\text{desorption + SCR + oxidation}}
$$

The volumetric NOx consumption rate involves only SCR, converted from surface to volumetric via $k_{s2v}$:

$$
\dot{[\text{NO}_x]}^{scr} = -k_{s2v}\,k_{scr}\,[\text{NH}_3]^{ads}\,[\text{NO}_x]^{in}
$$

**Observation.** These two rate equations are coupled through $[\text{NH}_3]^{ads}$ — the surface concentration drives NOx reduction but is itself unobservable from tailpipe measurements. The central challenge of this derivation is to eliminate this hidden state and express the NOx dynamics purely in terms of measurable signals.

---

## 2. Discretization

The reactor is modelled as a discrete plug-flow system. Within each sampling period $t_s$, there are $n$ residence times. At each residence time, a fresh parcel of exhaust gas enters and reacts with the catalyst surface. We now discretize each governing equation by tracking what happens across these $n$ sub-steps.

### 2a. Surface NH$_3$ dynamics

At each sub-step $i$, molar conservation on the catalyst surface gives:

$$
[\text{NH}_3]^{ads}(k{+}(i{+}1)\tau) = [\text{NH}_3]^{ads}(k{+}i\tau) + \tau\,\dot{[\text{NH}_3]}^{ads}(k{+}i\tau)
$$

where the rate is treated as constant within each residence time. The gaseous inputs ($[\text{NH}_3]^{in}$, $[\text{NO}_x]^{in}$) are held constant (ZOH) across the entire sample, so only the surface concentration changes sub-step to sub-step. Summing telescopically from $i = 0$ to $n{-}1$ (all intermediate terms cancel):

$$
[\text{NH}_3]^{ads}(k{+}1) - [\text{NH}_3]^{ads}(k) = \tau \sum_{i=0}^{n-1} \bigl(\Gamma\,\gamma_{ads}(k) - [\text{NH}_3]^{ads}(k{+}i\tau)\,\gamma_{des}(k)\bigr)
$$

where $\gamma_{ads} = k_{ads}\,[\text{NH}_3]^{in}$ and $\gamma_{des} = k_{ads}\,[\text{NH}_3]^{in} + k_{des} + k_{scr}\,[\text{NO}_x]^{in} + k_{oxi}$.

The sum over the unknown sub-step values $[\text{NH}_3]^{ads}(k{+}i\tau)$ defines the **average surface concentration** within the sample:

$$
\sigma(k) \;=\; \frac{1}{n}\sum_{i=0}^{n-1}[\text{NH}_3]^{ads}(k + i\tau)
$$

This $\sigma(k)$ is unobservable — we cannot measure the surface concentration at each sub-step. **Approximation (ZOH):** treat $\sigma(k) \approx [\text{NH}_3]^{ads}(k)$, i.e., the surface concentration at the start of the sample dominates the average. Using $n\tau = t_s$:

$$
\sigma(k{+}1) = \sigma(k) + t_s\,k_{ads}\,[\text{NH}_3]^{in}(k)\bigl(\Gamma - \sigma(k)\bigr) - t_s\,k_{od}\,\sigma(k) - t_s\,k_{scr}\,u_1(k)\,\sigma(k) \tag{S}
$$

### 2b. NOx output

Unlike adsorbed NH$_3$, NOx has no storage on the surface — the outlet concentration at each residence time depends only on the inlet concentration and the surface state during that sub-step:

$$
[\text{NO}_x]^{out}(k{+}(i{+}1)\tau) = [\text{NO}_x]^{in}(k{+}i\tau) - \tau\,k_{s2v}\,k_{scr}\,[\text{NH}_3]^{ads}(k{+}i\tau)\,[\text{NO}_x]^{in}(k{+}i\tau)
$$

There is no telescopic accumulation — the measured outlet NOx corresponds to the last residence time ($i = n{-}1$). Applying ZOH for inlet NOx ($u_1(k)$) and the same surface approximation ($\sigma(k)$):

$$
x_1(k{+}1) = u_1(k)\bigl(1 - \tau\,k_{s2v}\,k_{scr}\,\sigma(k)\bigr) \tag{N}
$$

### 2c. Urea dosing

The injected urea solution (AdBlue) undergoes evaporation, thermal decomposition, and isocyanic acid hydrolysis. These three steps aggregate to:

$$
\text{NH}_2\text{-CO-NH}_2 \;\longrightarrow\; 2\,\text{NH}_3 + \text{CO}_2 + x\,\text{H}_2\text{O}
$$

Since the urea concentration in AdBlue is fixed, the decomposition rate $r_u = k_u\,[\text{urea}]$ is a constant — it depends only on the volume of solution injected, not on temperature (the injector preheats the solution). The total moles of NH$_3$ produced in one sampling period from injected volume $t_s\,u_2(k)$:

$$
\{\text{NH}_3\}^{in}(k) = \underbrace{2\,r_u}_{\text{rate}} \times \underbrace{t_s\,u_2(k)}_{\text{volume injected}} \times \underbrace{t_s}_{\text{duration}}
$$

These moles are distributed evenly over the $n$ residence times. Dividing by the chamber volume $V$ and substituting $n = t_s/\tau$:

$$
[\text{NH}_3]^{in}(k) = \frac{\{\text{NH}_3\}^{in}(k)}{n\,V} = \frac{2\,r_u\,t_s\,\tau}{V}\;u_2(k)
$$

Substituting $\tau = V\rho_0/F(k)$ and defining $\nu_u = 2\,r_u\,t_s\,\rho_0$:

$$
[\text{NH}_3]^{in}(k) = \nu_u\,\frac{u_2(k)}{F(k)} \tag{U}
$$

The inlet NH$_3$ concentration is inversely proportional to flow — at lower flow, the same injection rate produces a higher concentration in the chamber. This reciprocal relationship breaks down at zero flow, which is physically meaningful (no exhaust to carry the reactant).

---

## 3. Eliminating the unobservable $\sigma$

Equations (S) and (N) form a coupled system, but $\sigma$ cannot be measured. The key idea is to express $\sigma$ in terms of an observable quantity derived from (N).

Define the **NOx reduction per sample**:

$$
\eta(k) = u_1(k{-}1) - x_1(k)
$$

From (N), $x_1(k{+}1) = u_1(k)(1 - \tau\,k_{s2v}\,k_{scr}\,\sigma(k))$, so $\eta(k{+}1) = \tau(k)\,k_{s2v}\,k_{scr}(k)\,u_1(k)\,\sigma(k)$. This gives:

$$
\sigma(k) = \frac{\eta(k{+}1)}{\tau(k)\,k_{s2v}\,k_{scr}(k)\,u_1(k)} \tag{E}
$$

Now substitute (E) into (S) to eliminate $\sigma$. Rewrite (S) compactly as:

$$
\sigma(k) = \sigma(k{-}1)\,\gamma_{proc}(k{-}1) + \Gamma\,t_s\,k_{ads}(k{-}1)\,[\text{NH}_3]^{in}(k{-}1)
$$

where $\gamma_{proc}(k) = 1 - t_s\bigl(k_{ads}(k)\,[\text{NH}_3]^{in}(k) + k_{scr}(k)\,u_1(k) + k_{od}(k)\bigr)$ collects all the terms that act on the current surface state. Applying (E) at times $k$ and $k{-}1$:

$$
\eta(k{+}1) = \eta(k)\;\frac{\tau(k)}{\tau(k{-}1)}\;\frac{u_1(k)}{u_1(k{-}1)}\;\frac{k_{scr}(k)}{k_{scr}(k{-}1)}\;\gamma_{proc}(k{-}1) \;+\; t_s\,k_{s2v}\,\Gamma\,\tau(k)\,u_1(k)\,[\text{NH}_3]^{in}(k{-}1)\,k_{scr}(k)\,k_{ads}(k{-}1)
$$

This recursion is entirely in terms of $\eta$ (observable) and the measurable inputs ($u_1, u_2, F, T$), but it still contains ratio prefactors and individual rate constants that need to be dealt with.

---

## 4. Simplification and parametrization

The recursion from Step 3 is exact but unwieldy. We now make a series of physically motivated simplifications that bring it into a form amenable to parameter estimation.

**Assumptions:**

| | Statement | Consequence |
|---|-----------|------------|
| A1 | $T(k) \approx T(k{-}1)$ (temperature changes slowly relative to sampling) | $k_{scr}(k)/k_{scr}(k{-}1) \approx 1$; also $k_{scr}(k)\,k_{ads}(k{-}1) \approx k_{scr/ads}(k{-}1)$, a combined rate constant with its own Arrhenius form |
| A2 | Density changes with temperature are negligible | $\tau(k) = \tau_0/F(k)$, so $\tau(k)/\tau(k{-}1) = F(k{-}1)/F(k)$ |
| A3 | Urea model (U) applies | $[\text{NH}_3]^{in}(k) = \nu_u\,u_2(k)/F(k)$ |
| A4 | Rate constants vary polynomially with $T$ over the operating range | Chebyshev parametrization (justified by Taylor expansion of Arrhenius) |
| A5 | $\Gamma$ is constant over the operating range (changes only with aging) | Absorbed into parameter vectors |

**Chebyshev temperature basis.** Rather than estimating Arrhenius parameters directly, the Arrhenius dependence is approximated by low-order polynomials in $T$. For numerical stability, the temperature range $[T_{min}, T_{max}]$ is mapped to $[-1,1]$ using Chebyshev polynomials ($T_0 = (T_{max}+T_{min})/2$, $T_r = (T_{max}-T_{min})/2$):

$$
\phi_1(k) = \begin{bmatrix} \tfrac{T-T_0}{T_r} & 1 \end{bmatrix}, \qquad
\phi_2(k) = \begin{bmatrix} 2\bigl(\tfrac{T-T_0}{T_r}\bigr)^2-1 & \tfrac{T-T_0}{T_r} & 1 \end{bmatrix}
$$

First-order $\phi_1$ is used for individual rate constants (nearly linear in $T$). Second-order $\phi_2$ is used for the product $k_{scr/ads}\,\Gamma$, which contains an inflection point.

**Flow-scaled NOx reduction.** After applying A1 and A2, the recursion still carries the ratio prefactor $\frac{u_1(k)}{u_1(k{-}1)}\frac{F(k{-}1)}{F(k)}$. To absorb this, define:

$$
\bar{\eta}_F(k) = \frac{F(k{-}1)}{u_1(k{-}1)}\;\eta(k)
$$

Multiplying the Step 3 recursion (with A1–A2 applied) by $F(k)/u_1(k)$ converts $\eta \to \bar{\eta}_F$ on both sides. The ratio prefactor becomes unity. Expanding $\gamma_{proc}$ and distributing, then substituting A3–A5 and the Chebyshev basis for each rate constant ($k_{ads} = \phi_1^T\theta_{ads}$, $k_{od} = \phi_1^T\theta_{od}$, $k_{scr} = \phi_1^T\theta_{scr}$, $k_{scr/ads} = \phi_2^T\theta_{scr/ads}$):

$$
\boxed{
\bar{\eta}_F(k{+}1) = \bar{\eta}_F(k)
+ u_{2F}(k{-}1)\,\phi_2(k{-}1)\,\theta_\Gamma
- u_{2F}(k{-}1)\,\bar{\eta}_F(k)\,\phi_1(k{-}1)\,\theta_{\eta_{ads}}
- \bar{\eta}_F(k)\,\phi_1(k{-}1)\,\theta_{\eta_{od}}
- u_1(k{-}1)\,\bar{\eta}_F(k)\,\phi_1(k{-}1)\,\theta_{\eta_{scr}}
}
\tag{D}
$$

where $u_{2F}(k) = u_2(k)/F(k)$, and the parameter vectors collect the physical constants:

| Parameter | Definition | Physical content |
|-----------|-----------|-----------------|
| $\theta_{\eta_{ads}}$ | $t_s\,\nu_u\,\theta_{ads}$ | site-blocking by adsorbed NH$_3$ |
| $\theta_{\eta_{od}}$ | $t_s\,\theta_{od}$ | oxidation + desorption loss |
| $\theta_{\eta_{scr}}$ | $t_s\,\theta_{scr}$ | SCR consumption of adsorbed NH$_3$ |
| $\theta_\Gamma$ | $t_s\,k_{s2v}\,\nu_u\,\Gamma\,\tau_0\,\theta_{scr/ads}$ | adsorption onto free sites ($\propto k_{scr/ads}\,\Gamma$) |

Each term in (D) has a clear physical interpretation. Reading right to left on the RHS: the current state persists ($\bar{\eta}_F(k)$), new NH$_3$ adsorbs onto free sites ($\theta_\Gamma$ term), and three loss channels drain the surface — site-blocking reduces the adsorption rate ($\theta_{\eta_{ads}}$), oxidation and desorption consume stored NH$_3$ ($\theta_{\eta_{od}}$), and SCR itself consumes stored NH$_3$ in proportion to inlet NOx ($\theta_{\eta_{scr}}$).

---

## 5. Catalyst saturation — the switched model

Equation (D) was derived assuming the catalyst operates in its unsaturated regime ($0 \le \sigma \le \Gamma$). In practice, the surface can saturate ($\sigma = \Gamma$, all sites occupied) or deplete ($\sigma = 0$, all sites empty). We need to impose these physical bounds.

**Translating the $\sigma$ constraint into $\bar{\eta}_F$ space.** From (E), $\sigma \propto \eta$, and $\bar{\eta}_F$ is a rescaled $\eta$. The constraint $0 \le \sigma \le \Gamma$ becomes bounds on $\bar{\eta}_F$:

- **Lower bound** ($\sigma = 0$): $\bar{\eta}_F = 0$ (no surface NH$_3$, no NOx reduction).
- **Upper bound** ($\sigma = \Gamma$): maximum possible $\bar{\eta}_F$ when every surface site is occupied.

**Saturated ceiling.** When $\sigma = \Gamma$, the NOx reduction is $\eta(k{+}1) = \tau(k)\,k_{s2v}\,k_{scr}(k)\,u_1(k)\,\Gamma$. Converting to $\bar{\eta}_F$ and parametrizing $k_{scr}$ with $\phi_2$ (quadratic, because $k_{scr} \times \Gamma$ has a non-monotone temperature dependence):

$$
\bar{\eta}_{F,sat}(k) = \phi_2^T(k)\,\theta_{sat}, \qquad \theta_{sat} = \Gamma\,\tau_0\,k_{s2v}\,\theta_{scr}
$$

$\theta_{sat}$ is a separate parameter vector (not shared with the unsaturated dynamics).

**Switched model.** Combining the unsaturated dynamics (D) with the saturation bounds:

$$
\boxed{
\bar{\eta}_F(k{+}1) = \max\!\bigl\{0,\;\min\!\bigl\{\text{(D)},\;\phi_2^T(k)\,\theta_{sat}\bigr\}\bigr\}
}
$$

| Condition | Regime | $\bar{\eta}_F(k{+}1)$ |
|-----------|--------|----------------------|
| $0 \le \text{(D)} \le \phi_2^T\theta_{sat}$ | Unsaturated | (D) |
| $\text{(D)} > \phi_2^T\theta_{sat}$ | Saturated ($\sigma = \Gamma$) | $\phi_2^T(k)\,\theta_{sat}$ |
| $\text{(D)} < 0$ | Depleted ($\sigma = 0$) | $0$ |

**Guard condition.** The catalyst remains unsaturated as long as the net adsorption onto free sites is non-negative, i.e., $\Gamma - \sigma \ge 0$. In (D), the adsorption-related terms are $u_{2F}\,\phi_2\,\theta_\Gamma - u_{2F}\,\bar{\eta}_F\,\phi_1\,\theta_{\eta_{ads}}$. Setting this $\ge 0$ and rearranging into a single inner product:

$$
\phi_g^T(k)\,\theta_{NO_x} \le 0, \qquad
\phi_g(k) = \begin{bmatrix} \bar{\eta}_F(k)\,\phi_1^T(k{-}1) \\ \mathbf{0} \\ \mathbf{0} \\ -\phi_2^T(k{-}1) \end{bmatrix}, \quad
\theta_{NO_x} = \begin{bmatrix} \theta_{\eta_{ads}} \\ \theta_{\eta_{od}} \\ \theta_{\eta_{scr}} \\ \theta_\Gamma \end{bmatrix}
$$

When $\phi_g^T\,\theta_{NO_x} > 0$, the system switches to saturated mode. This guard condition is *linear in the model parameters* — it can be evaluated from the current state and inputs without knowing $\sigma$. This structure — a threshold that depends on the state itself — makes the model a **Self-Excited Threshold Nonlinear ARX (SETNARX)** system.

**Tailpipe NOx** (converting back from $\bar{\eta}_F$):

$$
\eta(k{+}1) = \bar{\eta}_F(k{+}1)\;\frac{u_1(k)}{F(k)}, \qquad x_1(k{+}1) = u_1(k) - \eta(k{+}1)
$$

---

*Compiled for the SCR-dynamics JDSMC publication.*
