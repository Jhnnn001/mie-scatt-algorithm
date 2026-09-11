# plan.md — Analytic forward solver for a homogeneous spheroid (MATLAB)

Goal: a MATLAB function returning the complex electric field (Ex, Ey, Ez) at
arbitrary points for a plane wave scattered by a homogeneous prolate or oblate
spheroid, by the spheroidal-coordinate separation-of-variables solution
(Asano & Yamamoto 1975, Appl. Opt. 14, 29; Barton 2001, Appl. Opt. 40, 3598,
both in `ref/mie-scatt/`; spheroidal-function facts from DLMF Ch. 30 and
Flammer 1957). Intended use: ground truth for Born, Rytov and multi-slice BPM.

This document is self-contained: an implementer needs only this file, the two
papers, DLMF, and base MATLAB (no toolboxes). It lists every file, its purpose
and interface, the formulas, the numerical safeguards, and the tests that define
"done". It contains no full code. Revision history: v1 2026-09-12; v2 after
three review rounds (coordinate stability, radial-function stability,
convergence checks, validation scope); v3 implementation reconciliation
(arbitrarily normalized outgoing basis, fail-closed radial diagnostics,
rank-aware solves, runtime-validation scope); v4 post-implementation audit
(scaled oblate derivatives, overflow-safe geometry, complete test runner,
atomic stress reports); v5 final hardening (guaranteed-distance singular
nudges, fail-closed solution validation, and shared-unit and field-amplitude
invariance); v6 scalar-polarization scope characterization.

---

## 0. Decisions

| Topic | Decision |
|---|---|
| Shape | Spheroid centered at origin, symmetry axis = z. `a` = semi-axis along z, `b` = radius in the xy plane. `a>b` prolate, `a<b` oblate; `abs(a-b) <= 64*eps(max(a,b))` routes to the sphere (Mie fallback). |
| Units | Any length unit, shared by `lambda, a, b, X, Y, Z`. `lambda` = vacuum wavelength. |
| Media | Particle index `n_p` (complex allowed), background `n_m` (real). |
| Incidence | Arbitrary direction (`theta_inc`, `phi_inc`) and polarization `pol=[p_x,p_y]` (complex, not normalized). |
| Output | Default total field (outside: incident + scattered; inside: internal). `'scattered'` = total − incident everywhere (Born/Rytov definition); `'incident'`. |
| Accuracy | One knob `tol` (default 1e-6). Truncations, quadratures and ODE tolerances are chosen and verified automatically; failure to verify is an error by default. |
| Scope | A *runtime-validation envelope* (Section 6) with `info.validated` set only when the input is inside it **and** every runtime check passed. Outside it the solver attempts a solution subject to numerical ceilings, warns, and sets `validated=false`; exceptional or truncation-sensitive eigenproblems fail with an error. |
| Method | Spheroidal separation of variables with Galerkin (surface-projection) boundary matching (Barton). Incident field enters by numerical surface projection. |

---

## 1. Conventions (every file follows these)

1. Time dependence exp(−iωt). Plane wave E_inc(r) = E0 exp(i k_m k̂·r), phase 0 at origin.
2. k0 = 2π/λ, k_m = n_m k0, k_p = n_p k0. Maxwell in Asano's form (1): ∇×E = i k0 H, ∇×H = −i k0 n² E. Hence H = n k̂×E for a plane wave, and E = Σ(aN + bM) ⇒ H = −i n Σ(aM + bN) in the same region (Barton Eqs. 8–13).
3. Incidence: k̂ = (sinζ cosφ_i, sinζ sinφ_i, cosζ), ζ = `theta_inc`, φ_i = `phi_inc`. Rot = rotation taking ẑ to k̂ about ẑ×k̂ by ζ (identity at ζ=0). ê1 = Rot x̂, ê2 = Rot ŷ, E0 = p_x ê1 + p_y ê2. Closed form: ê1 = cosφ_i ê_TM − sinφ_i ê_TE, ê2 = sinφ_i ê_TM + cosφ_i ê_TE, ê_TE = (−sinφ_i, cosφ_i, 0), ê_TM = (cosζ cosφ_i, cosζ sinφ_i, −sinζ). Check: ζ=φ_i=0, [1,0] → x̂, [0,1] → ŷ.
4. Semifocal length l = √|a²−b²|, evaluated by `sph_semifocal.m` as `s·√[(|a−b|/s)(1+min(a,b)/s)]`, `s=max(a,b)`, so a shared rescaling cannot cause intermediate square overflow or underflow. All coordinate/function helpers work in units of l: c_ext = k_m l (real), c_int = k_p l (complex allowed), ξ0 = a/l for both types.
   - Prolate: x = √((ξ²−1)(1−η²)) cosφ, y = … sinφ, z = ξη, ξ ≥ 1. a/b = ξ0/√(ξ0²−1).
   - Oblate: x = √((ξ²+1)(1−η²)) cosφ, y = … sinφ, z = ξη, ξ ≥ 0, η signed with **sign(0) := +1** (upper sheet of the focal disk; either sheet gives the same field, see 2.1). b = l√(ξ0²+1).
   - (ê_η, ê_ξ, ê_φ) is right-handed in both systems.
5. Azimuthal index m ∈ [−M, M]; all spheroidal functions use μ = |m|; sign(m) appears only in e^{imφ} and in the i·sign(m) factors of φ-derivatives.
6. Scalar wavefunction Π_lm = S_lμ(c;η) R_lμ(c;ξ) e^{imφ}, l = μ..L. Angular functions normalized by ∫_{−1}^{1} S² dη = 1 (bilinear, also for complex c). Vector functions (lengths in l): M = ∇×(rΠ) = ∇Π × r, N = (1/c)∇×M.
7. Small geometric quantities α = ξ²−1 (prolate) or ξ²+1 (oblate), β = 1−η², D = ξ²∓η² are computed once, stably, in `sph_coords` and carried as fields; no helper recomputes them from ξ or η.

---

## 2. Mathematics

### 2.1 Coordinates (`sph_coords.m`) — cancellation-free recipes
Inputs x, y, z in units of l; ρ = hypot(x, y); φ = atan2(y, x).

**Prolate.** r_near = hypot(ρ, |z|−1), r_far = hypot(ρ, |z|+1) (distances to the focus on the same side of z=0 and the opposite one).
- ξ − 1 (stable): if |z| ≤ 1: (ρ²/2)[1/(r_far+1+|z|) + 1/(r_near+1−|z|)]; if |z| > 1: (|z|−1) + (ρ²/2)[1/(r_far+|z|+1) + 1/(r_near+|z|−1)]. ξ = 1 + (ξ−1).
- α = (ξ−1)(ξ+1); β = ρ²/α; D = α + β (= r_near r_far).
- η: if β ≥ 0.5 use η = sign(z)(r_far − r_near)/2 directly; else 1−|η| = β/(1+√(1−β)), η = sign(z)(1 − (1−|η|)). Both are accurate in their range; `beta` is the quantity used everywhere downstream.

**Oblate.** A = (ρ−1)(ρ+1) + z² (= ξ² − η²), D = hypot(A, 2z) (= ξ² + η²).
- If A ≥ 0, ξ = √((D+A)/2) and |η| = |z|/ξ. If A < 0, |η| = √((D−A)/2) and ξ = |z|/|η|. These product-identity forms avoid cancellation of the small coordinate.
- η = sign(z)|η| with sign(0) := +1; α = ξ² + 1; β = ρ²/α.
- Continuity across the disk (ξ=0, ρ<1) holds because R^(1)_lμ(ξ) has parity (−1)^{l−μ} and S_lμ(η) the same parity; either sheet gives identical E.

**Singular set.** Prolate foci (ρ=0, z=±1) and oblate focal ring (ρ=1, z=0) are singular points of the coordinate map (D→0). A point inside distance δ = 1e-8 (in units of l) is projected to distance δ without changing its nonsingular coordinate: prolate keeps z and sets ρ = √(δ²−(|z|−1)²), preserving φ when ρ>0 and choosing +x on the axis; oblate keeps z and sets ρ−1 to the same sign as its incoming value with magnitude √(δ²−z²), preserving φ and choosing the outward side when ρ=1 exactly. Thus a nudge never moves a point closer to the singular set. Documented accuracy at those points is ≈ c·δ. The axis, the prolate focal segment and the oblate disk are **not** nudged: all field formulas below are regular there. The exact axis is tested against nearby offsets, and nudged points against δ/2 on the same chosen side (T9).

**Returned struct `g`** (arrays over points): `eta, xi, phi, sigma, alpha, beta, D, nudged`, where `sigma=+1` for prolate and `sigma=-1` for oblate, `nudged` identifies projected singular-neighborhood inputs, `h_eta = √(D/β), h_xi = √(D/α), h_phi = √(αβ)`, and the vectors ê_η = (−η√(α/D) cosφ, −η√(α/D) sinφ, ξ√(β/D)), ê_ξ = (ξ√(β/D) cosφ, ξ√(β/D) sinφ, η√(α/D)), ê_φ = (−sinφ, cosφ, 0) (same form for both types with the type's α), plus the Euler-operator data of 2.5: A_ = ±ηβ (+ prolate, − oblate), A_' = ±(1−3η²), B = ξα, B' = 3ξ² ∓ 1, D_η = ∓2η, D_ξ = 2ξ. Downstream sign choices use `g.sigma`; the sign of A_ for oblate was derived from r·∇ = ξ̇∂_ξ + η̇∂_η under r→sr and is verified by T2.

### 2.2 Eigenvalues and Legendre coefficients (`sph_eigen.m`)
Angular ODE (Asano 13): d/dη[(1−η²)S'] + (λ ∓ c²η² − μ²/(1−η²))S = 0 (upper prolate).
Expand in orthonormal Ferrers functions P̄^μ_{μ+r}(η) (no Condon–Shortley phase), r ≡ l−μ (mod 2):
S_lμ = Σ'_r d̄_r P̄^μ_{μ+r}. Recurrence (Flammer 3.1.4 / DLMF 30.8) in the orthonormal basis is a **symmetric** tridiagonal eigenproblem:
- diagonal β_r = (μ+r)(μ+r+1) + [2(μ+r)(μ+r+1) − 2μ² − 1] c² / [(2μ+2r−1)(2μ+2r+3)]
- off-diagonal (r, r+2): e_r = c² √[(2μ+r+2)(2μ+r+1)(r+2)(r+1)] / [(2μ+2r+3)√((2μ+2r+1)(2μ+2r+5))]
- oblate: c² → −c².
Even-r and odd-r matrices separately; eigenvalues sorted by real part give l = μ, μ+2, … and μ+1, μ+3, …. Eigenvectors rescaled to Σ d̄_r² = 1 (bilinear), sign fixed by Re d̄_{l−μ} > 0.
- r_max: start at (L−μ) + ⌈|c|⌉ + 40 per parity class; recompute with 1.5 r_max and require ‖Δd̄‖ < 1e-13 for every retained l (independent of the L loop). Grow until satisfied.
- Complex c: labels may permute; the Galerkin solution is invariant to permutations *within the retained set*, so correctness rests on the cutoff checks of Section 3.3 (modal contribution tail and the (L+5, M+2) re-solve), not on labels. Additionally require the eigenvector matrix condition number < 1e8 and reject backward-stable eigenpairs that remain truncation-sensitive; exceptional or sensitive cases error even when their scalar parameters satisfy Section 6.
- c → 0 sanity: λ → l(l+1), d̄ → unit vector at r = l−μ.

### 2.3 Angular functions (`sph_angular.m`) — regular "reduced" form
Define reduced orthonormal functions Ũ^μ_ℓ(η) = P̄^μ_ℓ(η)/β^{μ/2} (polynomials in η).
- Seed Ũ^μ_μ = √((2μ+1)/2 · (2μ)!/(2^{2μ}(μ!)²)) (constant, via gammaln). Upward recurrence in ℓ at fixed μ (same as for P̄ because the β^{μ/2} factor is common): Ũ^μ_ℓ = a_ℓ[η Ũ^μ_{ℓ−1} − b_ℓ Ũ^μ_{ℓ−2}], a_ℓ = √((4ℓ²−1)/(ℓ²−μ²)), b_ℓ = √(((ℓ−1)²−μ²)/(4(ℓ−1)²−1)).
- Derivatives without division (DLMF 14.7 / Rodrigues): dŨ^μ_ℓ/dη = √((ℓ−μ)(ℓ+μ+1)) Ũ^{μ+1}_ℓ, d²Ũ^μ_ℓ/dη² = √((ℓ−μ)(ℓ+μ+1)(ℓ−μ−1)(ℓ+μ+2)) Ũ^{μ+2}_ℓ. These are applied **per Legendre term** ℓ = μ + r inside the sum, never per spheroidal mode l.
- S̃ = Σ_r d̄_r Ũ^μ_{μ+r}, S̃' and S̃'' with the per-term factors above (matrix products d̄ᵀ·Ũ).
- Outputs are the metric-weighted, everywhere-regular combinations (β from `g`, exponents ≥ 0 for μ ≥ 1; the S̃ terms vanish for μ = 0):
  - W0 = S = β^{μ/2} S̃
  - W1 = √β S' = β^{(μ+1)/2} S̃' − μη β^{(μ−1)/2} S̃
  - W2 = β^{3/2} S'' = β^{(μ+3)/2} S̃'' − 2μη β^{(μ+1)/2} S̃' + [μ(μ−2)η² − μβ] β^{(μ−1)/2} S̃
  - Wm = μ S/√β = μ β^{(μ−1)/2} S̃ (0 for μ = 0)
- Raw S' and S'' are never formed. Overflow note: Ũ^μ_ℓ(±1) ≈ √((2ℓ+1)(ℓ+μ)!/(2(ℓ−μ)!))/(2^μ μ!) stays < 1e150 for μ ≤ 300, ℓ ≤ μ+300; document μ ≤ 300 as a ceiling.

### 2.4 Radial functions (`sph_bessel.m`, `sph_radial.m`, `sph_radial_ode.m`)
First kind by the spherical-Bessel series (Asano 24–25, DLMF 30.11.3), valid for all ξ:
R^(1)_lμ = (α/ξ²)^{μ/2} Σ'_r i^{r+μ−l} d̄_r w̃_r j_{μ+r}(cξ) / Σ'_r d̄_r w̃_r, w̃_r = √[(2μ+2r+1)/(2μ+1) · C(r+2μ, r)] via gammaln. Terms with |d̄_r| < 1e-17 max|d̄| are dropped (safe for j only). Convergence is super-geometric.
- ξ ≥ 1 path (both types): j from `besselj` (√(π/2x) J_{n+1/2}), j' = j_{n−1} − (n+1)j/x, j'' = −(2/x)j' − (1 − n(n+1)/x²)j (spherical Bessel ODE, independent of the spheroidal ODE). Prefactor (α/ξ²)^{μ/2} evaluated with α from `g` and ξ² = α+1 (prolate) or α−1 (oblate), never from ξ directly.
- Oblate ξ < 1 path: with x = cξ define T_r(x) = j_{μ+r}(x)/x^μ (entire). Evaluate as exp(log|j| − μ log|x|)·sign(j) when j is finite and nonzero, 0 when j underflowed (then the term is provably negligible), and by the two-term series x^r/(2μ+2r+1)!!·(1 − x²/(2(2μ+2r+3))) when |x| < 1e-3. Derivatives: T_r' = (r/x)T_r − T_{r+1}, T_r'' = r(r−1)/x² T_r − (2r+1)/x T_{r+1} + T_{r+2} (the (r/x)T_r terms via the series form for tiny x). Then R^(1) = α^{μ/2} c^μ Σ(…)T_r(cξ)/Σ d̄w̃.
- Outputs are regular weighted combinations analogous to 2.3 (α from `g`; Ĝ = the r-sum without prefactor; prolate identity ξ² − α = 1):
  - V0 = R
  - V1 = √α R'. Prolate: = α^{(μ+1)/2} ξ^{−μ} Ĝ' + μ α^{(μ−1)/2} ξ^{−μ−1} Ĝ (from d ln(α/ξ²)^{μ/2}/dξ = μ(ξ²−α)/(αξ) = μ/(αξ), using ξ²−α = 1). Oblate: α ≥ 1 so R' has no singular factor; V1 = √α R' directly (small-ξ path via T_r').
  - V2 = α^{3/2} R'' expanded the same way; every exponent of α is ≥ 0 for μ ≥ 1, and the μ-terms vanish for μ = 0.
  - Va = R/√α = α^{(μ−1)/2} ξ^{−μ} Ĝ (μ ≥ 1; unused for μ = 0), Vm = μ Va
- **Outgoing exterior basis (real c_ext).** A canonical Neumann-series normalization is unnecessary for boundary matching and is numerically unreliable near the surface. For σ = +1 (prolate) or −1 (oblate), α = ξ² − σ, solve αR'' + 2ξR' + (c²ξ² − λ − σμ²/α)R = 0 for an outgoing basis R^(out), independently scaled per mode so that R^(out)(ξ0) = 1. This arbitrary scale is absorbed by the Galerkin coefficients; no canonical R^(2), splice, or Wronskian normalization is used.
  - Integrate the Riccati/log-amplitude state w = R'/R − ic + 1/ξ and q̃ inward with `ode113`, RelTol = 1e-12. For state ordering `[w; qtilde]`, AbsTol is `min(1,c)*1e-12` for each w component and `1e-12` for each q̃ component: w' = (λ−σc²)/α − 2σic/(αξ) + 2σ/(αξ²) + σμ²/α² − 2icw − w² − 2σw/(αξ), q̃' = w − 1/ξ. Reconstruct R(ξ) = exp[ic(ξ−ξ0) + q̃(ξ)−q̃(ξ0)] and R' = (ic − 1/ξ + w)R.
  - The outgoing seed is the adaptively truncated inverse-ξ expansion R = exp(icξ)ξ^−1Σ_n a_nξ^−n, with at most 80 coefficients. For c < 1, store b_n = a_n c^n and evaluate in q = 1/(cξ), avoiding coefficient overflow while retaining the derivative chain rule. Start with ξ_cmp = max(2ξ_switch, 3, c/3, √(max|λ|)/c), integrate from both ξ_cmp and ξ_start = 1.5ξ_cmp, and increase the starts by 1.5 (at most 12 attempts) until the normalized log-derivative and reconstructed-field differences are both < 1e-10.
  - Dense ODE output is used on [ξ0, ξ_start]; the same inverse-ξ expansion is used beyond ξ_start. Weighted outputs V0, V1, V2, Va and Vm are formed from R, R' and the radial ODE for R''.
  - Verification fails closed. At ODE-overlap and asymptotic points evaluate the normalized operator backward error | α(y' + y²) + 2ξy + c²ξ² − λ − σμ²/α | divided by | αy' | + | αy² | + | 2ξy | + | c²ξ² | + |λ| + |μ²/α|, where y = R'/R. The ODE and asymptotic errors and their numerical floors must all be ≤ 1e-9; otherwise `chk.ok=false`, and `spheroid_solve` rejects that μ basis.

### 2.5 Vector wavefunctions (`sph_vecwave.m`) — all terms regular
With e = e^{imφ}, s = sign(m), and W, V, g as above (± = + prolate, − oblate):
- g_η = W1 V0 e/√D, g_ξ = W0 V1 e/√D, g_φ = i s Wm Va e
- M_η = −i s Wm Va e ξ√(α/D); M_ξ = i s Wm Va e (±η√(β/D)); M_φ = e/D [ξ√α W1 V0 ∓ η√β W0 V1]
- Φ = ψ + r·∇ψ = e[W0V0 + (±η√β W1 V0 + ξ√α W0 V1)/D]
- N_η = (1/(c√D)) [W1V0 + (A_' W1V0 ± η W2V0 + ξ√α W1V1)/D − D_η(±ηβ W1V0 + ξ√(αβ) W0V1)/D²] e + c(±η√(β/D)) W0V0 e
- N_ξ = (1/(c√D)) [W0V1 + (B' W0V1 + ξ W0V2 ± η√β W1V1)/D − D_ξ(±η√(αβ) W1V0 + ξα W0V1)/D²] e + c ξ√(α/D) W0V0 e
- N_φ = (i s/c) [Wm Va + (±η W1 Vm + ξ Wm V1)/D] e
These come from M = ∇ψ×r and the identity ∇×∇×(rψ) = ∇(ψ + r·∇ψ) + c² rψ with r·∇ψ = (A_ψ_η + Bψ_ξ)/D; every metric factor has been absorbed so that only √α, √β, 1/√D, 1/D, η, ξ remain. The only remaining singularity is D → 0 (foci / focal ring, nudged). T2 checks ∇×M = cN, ∇×N = cM, ∇·M = ∇·N = 0 by finite differences, on and off the axis, on the focal segment and disk, for both types.

### 2.6 Field expansions
- Scattered (ξ ≥ ξ0, c_ext, R^(out)): E^(s) = Σ[a_lm N + b_lm M], H^(s) = −i n_m Σ[a M + b N].
- Internal (ξ < ξ0, c_int, R^(1)): E^(c) = Σ[c_lm N + d_lm M], H^(c) = −i n_p Σ[c M + d N].
- Incident: plane wave of Section 1; H^(i) = n_m k̂×E^(i).

### 2.7 Galerkin boundary equations (Barton 18–21, 26–29, one layer)
Surface nodes: Gauss–Legendre in θ ∈ [0, π] with η = cosθ and Jacobian sinθ (integrands are smooth in θ; in η they carry half-integer powers of β for odd μ, which would reduce Gauss–Legendre to algebraic convergence), Q_θ = L + max(c_ext, |c_int|) + 30 nodes; uniform φ nodes with FFT, Q_φ = 2(M + χ) + 32 with χ = |k_m b sinζ|.
For each m and test index l' = μ..L, with ⟨F⟩_{l'} := ∫ F(ξ0,η) S_{l'μ}(c_ext;η) dη and ⟨F^(i)⟩_{l'm} := (1/2π)∫∫ F^(i) S_{l'μ}(c_ext;η) e^{−imφ} dη dφ:
- E_η: Σ_l [a_l⟨N^(out)_η⟩ + b_l⟨M^(out)_η⟩ − c_l⟨N^(1)_η⟩ − d_l⟨M^(1)_η⟩] = −⟨E^(i)_η⟩; E_φ likewise.
- H_η: Σ_l [−i n_m(a_l⟨M^(out)_η⟩ + b_l⟨N^(out)_η⟩) + i n_p(c_l⟨M^(1)_η⟩ + d_l⟨N^(1)_η⟩)] = −⟨H^(i)_η⟩; H_φ likewise.
Tangential components of the incident field at the nodes come from E_inc·ê_η, E_inc·ê_φ (this multiplies the Jacobi–Anger spectrum J_m(χ(η)) by cos/sin(φ−φ_i), i.e. shifts it to J_{m±1}: the M margin below includes +1). Unknowns per m: 4(L−μ+1). Blocks assembled as (components at nodes)·diag(w_q sinθ_q)·S_extᵀ. m and −m are separate solves (ponytail: reflection symmetry could halve this).

### 2.8 Evaluation (`spheroid_eval.m`)
For `side='auto'` (default), classify directly in Cartesian coordinates: (X/b)² + (Y/b)² + (Z/a)² ≥ 1 − 64 eps is exterior, so the boundary tie is exterior without a cancellation-prone coordinate conversion. `'in'`/`'out'` force a side and are allowed only for |ξ−ξ0| ≤ 1e-8 (else error). Exterior: E^(i) + E^(s) using the outgoing ODE basis through ξ_start and its inverse-ξ continuation beyond; interior: E^(c). `'scattered'` returns the exterior expansion directly, without cancellation through `(E_inc + E_sca) - E_inc`, and returns E^(c) − E^(i) inside. `'incident'` is an analytic fast route in `spheroid_field` and does not invoke the scattering solve, while retaining forced-side validation. An empty `mu_data` sector is accepted only when every corresponding mode is analytically zero and all four coefficient vectors are finite, exactly zero, and have length L−μ+1; malformed sectors raise `spheroid_eval:InvalidSolution`. Convert each point chunk to `g`, then transform components to Cartesian with ê. Returns E and H.

### 2.9 Sphere (`mie_field.m`)
Axes differing by at most 64 floating-point spacings: Bohren & Huffman Mie (exp(−iωt), n'' ≥ 0 absorbing — the same conventions as this solver). A wider fixed near-sphere approximation is not used because it cannot honor arbitrary `tol` or the actual spheroid boundary-side contract. The implementation uses a_n, b_n, c_n, d_n; vector spherical harmonics with π_n, τ_n; and the rotated frame of Section 1. The automatic sphere-boundary tie treats r/a ≥ 1 − 64 eps as exterior. Its coefficient/sum convergence criterion is fixed at 1e-10; a stricter requested `tol` still returns the computed Mie field but sets `info.validated=false` and records the failed `mie_tolerance` check rather than claiming that target. Benchmark fixed now (T3): verified Wiscombe MVTstNew case 14, converted from its opposite absorption convention to m = 1.5 + i for exp(−iωt), at x = 1; expected Q_ext = 2.336321, Q_sca = 0.6634538, and Q_back = 0.573002538306. This verified fixture supersedes the unverified BHMIE Appendix-A digits previously recalled from memory. The test also asserts the analytic Rayleigh limit Q_sca = (8/3)x⁴|(n²−1)/(n²+2)|² at x = 0.01 to 1e-4 relative, and the optical theorem C_ext = (2π/k²)Σ(2n+1)Re(a_n+b_n).

---

## 3. Numerical safeguards

### 3.1 Start values
- x_ext = k_m max(a,b), x_int = |k_p| max(a,b), x = max(x_ext, x_int); L0 = ⌈x + 4.05 x^{1/3} + 10⌉.
- χ = |k_m b sinζ|; M0 = ⌈χ + 4χ^{1/3} + 10⌉ + 1 (the +1 for the J_{m±1} shift), M0 ≥ 1, M0 ≤ L0.
- Q_θ, Q_φ as in 2.7; r_max as in 2.2; ξ_switch, ξ_start as in 2.4.

### 3.2 Column equilibration
The outgoing basis has R^(out)(ξ0)=1, while internal regular functions and vector derivatives can still produce strongly unequal column scales. Divide every assembled column by its max-abs. At axial incidence (within floating-point angular resolution), construct only the |m|=1 sector; all omitted `mu_data` sectors are empty and their modal coefficients are exact zeros. Also preserve an exactly zero RHS as zero. Otherwise use `lsqminnorm` when rcond ≤ N eps for the N×N equilibrated matrix and `\` elsewhere, then undo the column scaling. Report minimum active rcond, numerical rank and solver method in `info`.

### 3.3 Convergence loop (driven by `tol`)
After each solve, all of the following must hold; otherwise refine the named quantity and re-solve (at most 5 rounds, then error unless `opts.on_fail='warn'`):
1. **Incident spectrum tail**: max over the last three shells |m| ∈ {M−2, M−1, M} of ‖RHS_m‖ ≤ tol·max_m‖RHS_m‖ (three shells because J_m has zeros). Else M ← M + ⌈0.25M⌉ + 2.
2. **Q_φ independence**: recompute the RHS with 2Q_φ; relative change ≤ tol. Else Q_φ ← 2Q_φ. (Checked separately from M.)
3. **Modal contribution tail** (normalization-free): for each m, the field contributed on the surface check grid by the top three retained l (|a N + b M| and the H analog, and likewise for c,d on the inside) ≤ tol·|E0|. Else L ← ⌈1.25 L⌉ + 5. Coefficient magnitudes are *not* used (they depend on scaling and normalization).
4. **Surface residual**: on an independent check grid (Chebyshev η, L+10 points × uniform φ, 2(2M+5) points, offset half a step from the FFT nodes) the jumps of tangential E, tangential H and n²E_ξ across ξ0, relative to |E0| and n_m|E0|, are ≤ tol. Then double the grid in both directions: the residual must not grow by more than 10% (Fourier tail stable). Else L, M as above and Q_θ ← ⌈1.5 Q_θ⌉.
5. **Q_θ independence**: re-solve with 1.5 Q_θ; compare the reconstructed E, H at 50 fixed random points (inside and outside) and C_sca: change ≤ tol. Coefficients are not compared.
6. **Cutoff stability**: re-solve at (L+5, M+2): E, H at the same 50 points change ≤ tol. This is what protects against mode swaps at the cutoff for complex c.
7. **Special-function checks** (2.2 r_max and eigenvector conditioning, 2.4 independent-start agreement and normalized ODE/asymptotic backward errors) all passed.
`info` reports L, M, Q_θ, Q_φ, r_max, ξ_switch, ξ_start, residual, the minimum rcond over active sectors, every check's value, and `validated`. Per-μ diagnostic arrays use NaN for analytically omitted axial sectors rather than treating them as solved values.

### 3.4 `info.validated`
For the spheroidal route, true iff the input lies in the runtime-validation envelope (Section 6) **and** all checks in 3.3 passed. Inside-envelope alone is not sufficient; envelope membership without passing checks yields false. The sphere route instead applies the fixed Mie tolerance rule in 2.9.

### 3.5 Overflow/underflow
Orthonormal reduced Legendre (no (2μ)! growth), gammaln for w̃_r and seeds, log-form T_r, dropping negligible d̄_r in the j-series only, chunked Bessel evaluation. Cartesian radii use nested `hypot`; surface normals are formed from dimensionless axis ratios; convergence cross sections square the far-field amplitude only after division by |E0| and length normalization. The same ordering avoids intermediate square underflow or overflow at the tested shared length factors 1e-200 and 1e200 and the tested spheroidal polarization amplitudes 1e-200 and 1e200. Ceilings: μ ≤ 300, L ≤ 400 (Ũ and C(ℓ+μ, 2μ) stay below 1e300).

### 3.6 Cost
Assembly + solves ≈ Σ_m O((4(L−μ+1))³); evaluation scales with points × Σ_m L and the angular/radial work. No wall-clock estimate is claimed before the coverage campaign is completed. ponytail: no `parfor` initially; evaluation chunks are the obvious place to add it.

---

## 4. Files and interfaces (`spheroid-analytic-forward/`)

User-unit interfaces: `spheroid_field`, `spheroid_solve`, `spheroid_eval`, `incident_plane_wave`, `mie_field`. Units of l inside `sph_*`.

| File | Purpose | Interface |
|---|---|---|
| `spheroid_field.m` | Entry point: validate, choose prolate/oblate/Mie, solve once, evaluate, reshape. | `[Ex,Ey,Ez,info] = spheroid_field(lambda,n_p,n_m,a,b,theta_inc,phi_inc,pol,X,Y,Z,opts)`; `opts`: `tol` (1e-6), `field`, `side`, `Lmax`, `Mmax`, `on_fail` ('error'), `verbose`. `info`: everything listed in 3.3 plus `sol`, `type`, `xi0`, `c_ext`, `c_int`, `time`. |
| `spheroid_solve.m` | Geometry, start values, quadrature, assembly (2.7), RHS projection, equilibrated solve, the loop of 3.3, `validated`. | `sol = spheroid_solve(lambda,n_p,n_m,a,b,theta_inc,phi_inc,pol,opts)`. |
| `spheroid_eval.m` | Section 2.8. | `[E,H] = spheroid_eval(sol,X,Y,Z,field,side)`, E,H of size [numel(X),3]. |
| `sph_semifocal.m` | Overflow-safe semifocal length from two semi-axes. | `l = sph_semifocal(a,b)`. |
| `sph_coords.m` | Section 2.1 including nudges and all returned fields. | `g = sph_coords(type,x,y,z)`. |
| `sph_eigen.m` | Section 2.2 with r_max convergence and cond check. | `eg = sph_eigen(type,mu,c,L)` → `lambda`, `d` (r × l), `r`, `rmax`, `cond`. |
| `sph_angular.m` | Section 2.3, outputs W0, W1, W2, Wm (l × points). | `W = sph_angular(eg,g)`. |
| `sph_bessel.m` | Spherical j, y and derivatives for order vectors × argument vectors (chunked). | `[z,dz,d2z] = sph_bessel(kind,nu,x)`. |
| `sph_bessel_T.m` | Log-safe T_r(x) = j_{μ+r}(x)/x^μ and derivatives (2.4 oblate small-ξ path), with the small-x series. | `[T,dT,d2T] = sph_bessel_T(mu,r,x)`. |
| `sph_radial.m` | Section 2.4 regular first-kind series, with a regular-IVP fallback when its normalization sum is ill-conditioned; also exposes diagnostic second/third-kind series for ξ>1. Returns V0, V1, V2, Va, Vm and `ok`. | `[V,ok] = sph_radial(eg,kind,g,c)`. |
| `sph_radial_ode.m` | Section 2.4 arbitrarily normalized outgoing basis: inverse-ξ seed, transformed inward ODE, dense output, independent-start and fail-closed backward-error checks. | `[rod,chk] = sph_radial_ode(eg,c,xi0,xi_switch)`. |
| `sph_radial_ode_eval.m` | Evaluate the stored outgoing ODE solution or its asymptotic continuation. | `V = sph_radial_ode_eval(rod,g)`. |
| `sph_vecwave.m` | Section 2.5. | `[M,N] = sph_vecwave(W,V,g,m,c)` → structs with `eta, xi, phi` (l × points), e^{imφ} excluded. |
| `incident_plane_wave.m` | Section 1: E_inc, H_inc, k̂, E0. | `[E,H,khat,E0] = incident_plane_wave(k_m,n_m,theta_inc,phi_inc,pol,X,Y,Z)`. |
| `mie_field.m` | Section 2.9. | `[E,H] = mie_field(lambda,n_p,n_m,a,theta_inc,phi_inc,pol,X,Y,Z,field,side)`. |
| `mie_efficiencies.m` | Mie efficiencies and optional coefficient vectors for the sphere route. | `[Qext,Qsca,Qback,coeff] = mie_efficiencies(x,n_rel)`. |
| `gauss_legendre.m` | Golub–Welsch nodes/weights. | `[x,w] = gauss_legendre(Q)`. |
| `example_xz_slice.m` | Demo: |E| on the x–z plane for the Barton prolate case. | script |
| `tests/test_gauss_legendre.m` | Quadrature exactness and validation. | function with asserts |
| `tests/test_incident_plane_wave.m` | Plane-wave conventions, Maxwell relation and validation. | function with asserts |
| `tests/test_sph_coords.m` | Coordinate identities, singular nudges and semifocal scaling. | function with asserts |
| `tests/test_solution_validation.m` | Fail-closed evaluator validation for malformed active sectors. | function with asserts |
| `tests/test_unit_scaling.m` | Mie and spheroidal invariance under shared length factors 1e-200 and 1e200; spheroidal invariance also under field-amplitude factors 1e-200 and 1e200. | function with asserts |
| `tests/test_special_functions.m` | T1 | function with asserts |
| `tests/test_coords_vecwave.m` | T2, T9 | function with asserts |
| `tests/test_sphere_mie.m` | T3 (Mie benchmark, near-sphere) | function with asserts |
| `tests/test_scalar_approximation.m` | T10 fixed-polarization scalar scope | function with asserts |
| `tests/test_physics.m` | T5–T8 | function with asserts |
| `tests/test_barton2001.m` | T4 | function with asserts |
| `tests/stress_envelope.m` | Section 6 coverage report (not pass/fail) | script |
| `tests/run_all.m` | Runs all eleven test functions, prints pass/fail and timings. | script |

Dependencies: base MATLAB R2018b+ (`besselj`, `bessely`, `gammaln`, `eig`, `fft`, `ode113`, `deval`, `lsqminnorm`, and the `max(...,'all')` syntax). No toolboxes.

---

## 5. Tests (definition of done)

T1 `test_special_functions.m`
- c = 1e-6: λ − l(l+1) < 1e-8; d̄ unit vectors; R^(1) ≈ j_l(cξ), R^(2) ≈ y_l(cξ) at cξ = 3 to 1e-8.
- Orthonormality ∫S S dη = δ to 1e-10, (prolate, oblate) × c ∈ {0.5, 5, 50, 200} × μ ∈ {0, 3, 40, 200}.
- Angular ODE residual using W-based S, S', S'' and regular-radial ODE residual using independently formed derivatives are < 1e-9. The outgoing regression cases `(type,c,μ,L,ξ0)` are `(prolate,50,0,2,1.03)`, `(prolate,200,0,2,1.03)`, `(prolate,50,40,42,1.03)`, `(prolate,200,200,200,1.03)`, `(oblate,50,0,2,0.577)`, `(oblate,200,0,2,0.26)`, `(oblate,50,40,42,0.26)`, and `(oblate,200,200,200,0.577)`. Each requires R(ξ0)=1, independent-start log-derivative and field differences < 1e-10, normalized ODE/asymptotic backward errors and their floors ≤ 1e-9, finite far evaluation, and an independent finite-difference radial residual < 2e-7.
- The oblate inner path checks V0, V1 and V2 against the scaled T_r power series just above |cξ|=1e-3, where the unscaled derivative differences lose digits; r_max independence (2.2) and eigenvector cond are also reported.

T2 `test_coords_vecwave.m`
- Finite-difference checks ∇×M = cN, ∇×N = cM, ∇·M = ∇·N = 0 (step 1e-5, relative error < 1e-6) in at least 30 type/mode/basis cases including: exact z-axis (inside and outside), prolate focal segment, oblate disk (both signs of z at |z| = 1e-9), generic points; both types, m ∈ {−3, 0, 2, 7}, regular/outgoing series, and the production outgoing ODE basis.
- T9 (singularities): field at exact axis points vs offsets ρ = 1e-6, 1e-7 differs ≤ 3c·ρ|E0|; nudged focus/ring evaluations with δ and δ/2 agree ≤ 1e-7|E0|; continuity across the oblate disk (z = ±1e-9) ≤ 1e-8|E0|. Separately, `test_sph_coords.m` checks α and β at r = 1e7 l against their direct spherical-limit formulas to 1e-12 relative.

T3 `test_sphere_mie.m`
- Verified Wiscombe MVTstNew case 14 of 2.9, Rayleigh limit, optical theorem, tangential continuity at r = a to 1e-8, and the fixed Mie tolerance report.
- Near-sphere: a/b = 1 ± 1e-4, k_m b = 10, n_p/n_m = 1.2, 41×41 grid in the plane y = 0.3b; cases (ζ=0, pol=[1,0]) and (ζ=40°, φ_i=70°, pol=[1, 0.5i]); max|E_spheroid − E_Mie|/|E0| < 1e-2. Axes differing by at most 64 floating-point spacings route to Mie exactly.

Auxiliary regressions: `test_sph_coords.m` verifies semifocal lengths at scales 1e-200, 1 and 1e200 and both-sided guaranteed-distance nudges; `test_solution_validation.m` corrupts active-sector metadata and stored ODE payloads and requires `spheroid_eval:InvalidSolution`; `test_unit_scaling.m` compares Mie and spheroidal E,H under shared length factors 1e-200 and 1e200, then separately compares spheroidal E,H under polarization-amplitude factors 1e-200 and 1e200.

T4 `test_barton2001.m` (Barton 2001 Figs. 7, 13; n_m = 1, ζ = 30°, φ_i = 0, pol = [0,1] i.e. E ⊥ x–z plane; l = 1, λ = 2π/h_ext)
- Prolate 2:1: ξ0 = 1.154701, h_ext = 27.494593, a = ξ0, b = √(ξ0²−1). Oblate 2:1: ξ0 = 0.577350, h_ext = 21.822473, a = ξ0, b = √(ξ0²+1) (exercises the ODE path, ξ0 < 1).
- Scattered field on the x–z plane at r = 1e7 l, θ step 0.25°, S_r = (r²/(π l²)) n_m |E_s|² (Barton Eq. 50 with H = n k̂×E in the far zone). Compare the direct sampled maximum; a continuous-angle parabolic refinement is a different observable and is not used.
- Expected sampled maxima: prolate n_p = 1.33 → 12.18744, 1.50 → 14.48922; oblate 1.33 → 60.93952, 1.50 → 51.01663. Relative error < 1e-4. `info.validated` must be true.

T5–T8 `test_physics.m`
- T5 zero contrast (n_p = n_m, prolate and oblate, ζ = 35°): the analytic identity route gives |E_sca| ≤ 1e-8|E0| at 100 random points and internal = incident to 1e-8. It intentionally does not assemble the singular zero-contrast Galerkin system.
- T6 Rayleigh limit (x = 0.02, n_p/n_m = 1.5, both types): the quasistatic leading field is E_QS,j = E0,j/(1 + L_j(ε_r − 1)), ε_r = (n_p/n_m)². Check E_QS at the center and against the inversion-even average [E(r)+E(−r)]/2 at three interior point pairs; the average removes the O(kr) plane-wave phase gradient, leaving an O(x²) comparison with tolerance 1e-3. Prolate (e² = 1 − b²/a²): L_z = ((1−e²)/e²)[(1/(2e)) ln((1+e)/(1−e)) − 1], L_x = L_y = (1−L_z)/2. Oblate (e² = 1 − a²/b²), g = √((1−e²)/e²): L_x = L_y = (g/(2e²))[π/2 − atan g] − g²/2, L_z = 1 − 2L_x (Bohren & Huffman 5.33–5.34). E ∥ x via (ζ=0, [1,0]); E ∥ z via (ζ=90°, φ_i=0, [1,0]) (gives −ẑ).
- T7 energy: first calibrate the far-field formulas on the Mie sphere: C_sca = (1/|E0|²)∮ r²|E_s|² dΩ at r = 1e7 l (Gauss θ × uniform φ) and the optical theorem C_ext = (4π/k_m) Im[f(k̂)·ê0*]/|E0|² with E_s → f e^{ikr}/r must match the Mie series values (fixes signs). Then lossless spheroid (n_p = 1.5, x = 30, both types, ζ = 30°): C_ext = C_sca to 1e-6. Lossy spheroid (n_p = 1.4 + 0.02i): C_ext = C_sca + C_abs with C_abs computed two independent ways, C_abs^P = −(1/(n_m|E0|²))∮ Re(E×H*)·n̂ dA on a sphere enclosing the spheroid, and C_abs^V = (k0 Im(n_p²)/(n_m|E0|²)) ∫_V |E_in|² dV (3-D Gauss quadrature in η, ξ, φ with the volume element h_η h_ξ h_φ); all three relations to 1e-5.
- T8 boundary traces: with `side='in'` and `'out'` on 2(L+10)×2(2M+5) surface points, tangential E, tangential H and n²E_ξ continuity ≤ tol; report the values.

T10 `test_scalar_approximation.m`
- This characterizes only the fixed-polarization reduction, not the amplitude or phase accuracy of a particular scalar wave solver. For unit incident polarization ê0, define E_co = (E_s·ê0*)ê0 and leakage ε_perp = ||E_s−E_co||_F/||E_s||_F on generic off-symmetry far-field directions.
- The cell-sized cases use vacuum wavelength 0.532 µm, water n_m = 1.335381534 at 20 °C, a homogeneous nonabsorbing cell n_p = 1.37, and semi-axes (a,b) = (5,2.5) µm prolate or (2.5,5) µm oblate. Thus x_ext = 78.86. Global x and y polarization propagate along +z; global z polarization propagates along +x because a plane wave requires E0·k̂ = 0.
- At r = 1e4 max(a,b), directions proportional to k̂ + u ê0 + v(k̂×ê0) use nonzero staggered u,v grids. These are gnomonic square angular patches, not solid-angle-weighted cones. The narrow patch |u|,|v| ≤ 0.08 requires ε_perp < 0.1; the wide patch |u|,|v| ≤ 1 requires ε_perp > 0.1. Thus the test accepts a narrow-angle fixed-polarization reduction and rejects a global fixed-polarization claim.

Coverage `stress_envelope.m` (separate resumable report, not pass/fail): corners of Section 6, 20 Latin-hypercube interior samples, and a resonance sweep (n_p/n_m = 1.5, a/b = 1.2, x_ext from 30 to 35 in steps of 0.01); records residual, all 3.3 checks, time, `validated`. The MAT report is authoritative and promoted before its derived CSV; resume rejects duplicate, foreign or count-inconsistent rows. A stub-solver smoke verifies all 537 case definitions and this report/resume behavior. Regression coverage of a continuous domain, not a proof; this plan does not claim that the full numerical campaign has completed or passed.

---

## 6. Runtime-validation envelope

| Quantity | Scalar eligibility bound |
|---|---|
| x_ext = k_m max(a,b) | 0.02–250 |
| x_int = |k_p| max(a,b) | ≤ 250 |
| a/b | 1/4 … 4 (ξ0 ≥ 1.03 prolate; ξ0 ≥ 1/√15 ≈ 0.258 oblate) |
| Re(n_p/n_m) | 0.5 … 2.0 |
| Im(n_p)/Re(n_p) | ≤ 0.1 |
| M, L | ≤ 300, ≤ 400 |
These bounds determine eligibility for `info.validated`; they do not certify every point in the continuous rectangle. Exceptional or truncation-sensitive eigenproblems are rejected, including any such complex-oblate case within the scalar bounds. Outside the bounds the solver attempts a solution subject to numerical ceilings, warns and sets `validated=false`; inside, `validated` still requires every runtime check (3.4). The separate `stress_envelope` campaign supplies empirical coverage only after its report has actually completed.

---

## 7. Implementation order (each step ends with its verification)

1. `gauss_legendre`, `sph_coords`, `incident_plane_wave`, `mie_field`, `mie_efficiencies` → T3 Mie parts (verified Wiscombe MVTstNew case 14), finite-difference ∇×E_inc = i k0 H_inc, T9 coordinate checks (α, β, D against direct formulas away from singular sets).
2. `sph_eigen`, `sph_bessel`, `sph_angular`, `sph_radial`, `sph_radial_ode` → T1.
3. `sph_vecwave` → T2 (including axis, focal segment, disk).
4. `spheroid_solve`, `spheroid_eval`, `spheroid_field` with the full 3.3 loop → near-sphere part of T3, T5, T8.
5. T4, T6, T7, `example_xz_slice`, `stress_envelope`, `run_all` → the original ten-suite pass/fail run, no-JVM example, and 537-case stress-report smoke are verified. T10 is verified separately. The full numerical stress campaign remains a separate unexecuted coverage study and is not claimed below.

Deliberate simplifications (each marked `ponytail:` in code): −m solved separately; no parallel evaluation; plane-wave illumination only (a focused beam changes only the incident projection and `incident_plane_wave`).

---

## 8. Verification result (2026-09-12)

- MATLAB R2024a, launched with `-nojvm`: the ten pre-T10 functions passed at commit `569307b`; the authoritative log is `../tmp/spheroid-analytic-forward/green/run_all_final10.log`. `test_physics` took 7520.810 s and `test_barton2001` took 1594.542 s. T10 is verified separately after each parameter revision. The complete eleven-suite runner has not been rerun because its two slowest unchanged suites require about 2.5 hours combined.
- The largest lossless energy mismatch was 9.89e-8. The lossy surface-flux and volume-loss absorption values met the 1e-5 criterion; refining the volume quadrature changed its result by 5.23e-16 (prolate) and 3.18e-15 (oblate). The largest T8 boundary component was 1.86e-7.
- The four Barton sampled-maximum relative errors were 7.92e-8, 5.66e-7, 3.02e-7 and 4.25e-7, all below 1e-4.
- `example_xz_slice.m` evaluated 90,601 finite field magnitudes under `-nojvm` and skipped plotting as designed; the range was 0.00432103841629–5.97273629666. `checkcode` reported zero issues across all 30 MATLAB files.
- The stress-report stub smoke verified all 537 case definitions and atomic/resumable report behavior. The full 537-case numerical campaign has not been run, so it supplies no additional envelope-certification claim.
