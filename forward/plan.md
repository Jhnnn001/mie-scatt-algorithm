# plan.md — Born and Rytov forward models for a homogeneous spheroid (MATLAB)

Goal: a MATLAB function returning the first-order Born and Rytov scattered
fields u_B, u_R (scalar) at arbitrary points for a plane wave incident on a
homogeneous prolate or oblate spheroid, and a comparison script that measures
both approximations against the vector solution of
`../spheroid-analytic-forward/` (the ground truth).

References (all in `forward/ref/`):
- Kak & Slaney, *Principles of Computerized Tomographic Imaging*, Ch. 6
  (`SIAMB0000017_chapter-10_...pdf`): §6.1.2 Eqs. 28–37 (inhomogeneous
  Helmholtz, Green's function, integral equation), §6.2.1 Eqs. 38–51 (Born,
  validity), §6.2.2 Eqs. 52–86 (Rytov, validity, u_B ↔ φ_s).
- Müller, Schürmann, Guck 2016, *The Theory of Diffraction Tomography*,
  arXiv:1507.00466v3 (`1507.00466v3.pdf`): §2.2 Eqs. 2.8–2.16, §3.1 Eqs.
  3.1–3.15 (Born), §3.2 Eqs. 3.16–3.52 (Rytov), §5.1 (3-D Fourier
  diffraction theorem, used only as background).
- Sibling `../spheroid-analytic-forward/plan.md` for geometry, incidence and
  the exp(−iωt) convention, which this folder shares exactly.

This document is self-contained: an implementer needs only this file, the two
references and base MATLAB. It lists every file, interface, formula, safeguard
and the tests that define "done". It contains no full code.
Revision: v1 2026-09-12; v2 2026-09-12 after three review rounds (phase
unwrapping criteria, far-field kernel and accuracy statement, `expm1`,
non-finite Rytov output, N_max rule and surface tests, vector-Rytov wording,
ground-truth status); v3 2026-09-12 after implementation verification
(independent Born/Rytov convergence, final N_max cap pair, physical-clearance
envelope, reference-quadrature order, ground-truth routing, and bounded chunk
memory).

---

## 0. Decisions

| Topic | Decision |
|---|---|
| Model | Scalar Helmholtz (Kak 28, Müller 2.8). The scalar u stands for the co-polarized component of E (Kak p. 205: "u may be set equal to the complex amplitude of the electric field along its polarization"). No vector Born/Rytov (Section 8). |
| Object | Homogeneous spheroid only, same parametrization as the sibling: centred at the origin, symmetry axis z, `a` = semi-axis along z, `b` = radius in xy, `a==b` sphere. `n_p` complex allowed, `n_m` real. No voxel grids. |
| Incidence | Plane wave, direction (`theta_inc`, `phi_inc`) exactly as sibling §1.3; unit amplitude u0 = exp(i k_m k̂·r). No polarization argument (scalar). |
| Output | u_B, u_R and u0 at arbitrary points (X, Y, Z), inside or outside the spheroid. Scattered fields only; total fields are u0+u_B and u0+u_R (Kak 38, 81). |
| Numerics | The Born volume integral over the spheroid is reduced *exactly* to a surface integral over the spheroid (Green's second identity, §2.1) plus a jump term. 2-D Gauss–Legendre(θ) × trapezoid(φ) quadrature, spectral convergence, no volume singularity for interior points. Convergence verified automatically per point (`tol`). |
| Rytov | u_R = u0 · expm1(u_B/u0) (Kak 72 + 84, Müller 3.37–3.38). Forward direction needs no phase unwrapping. |
| Dependencies | Base MATLAB. `../spheroid-analytic-forward/` on the path for `gauss_legendre.m`, `incident_plane_wave.m` (k̂, E0 in the comparison script) and, for Section 6 only, `mie_field.m` / `spheroid_field.m` (ground truth; see the status note in Section 6). |
| Accuracy | One knob `tol` (default 1e-8, relative per point, applied independently to the Born and Rytov fields with the common phase e^{ik\|r\|} removed — 2.3, 3.2). Points that fail to converge within `N_max`, or whose output is not finite, raise an error by default; `opts.on_fail='warn'` returns them with a mask and `validated=false`. |

---

## 1. Conventions (every file follows these)

1. Time dependence exp(−iωt). k0 = 2π/λ, k := k_m = n_m k0 (real). Plane wave u0(r) = exp(i k k̂·r), k̂ = (sinζ cosφ_i, sinζ sinφ_i, cosζ), ζ = `theta_inc`, φ_i = `phi_inc` (sibling §1.3; `incident_plane_wave` returns the same k̂). Unit amplitude ⇔ |E0| = 1.
2. Scattering potential (Müller 2.9; Kak 29 with n → n_p/n_m): f = k0² (n_p² − n_m²) = k² [(n_p/n_m)² − 1], constant inside V, zero outside. Complex for absorbing n_p.
3. Green's function G(R) = exp(ikR)/(4πR), (∇² + k²)G = −δ (Kak 32–33, Müller 2.10–2.11). Outgoing for exp(−iωt), consistent with the sibling's R^(3) choice.
4. First Born (Kak 40, Müller 3.6): u_B(r) = f ∫_V G(|r − r'|) u0(r') d³r'.
   First Rytov (Kak 71–72, 84; Müller 3.36–3.38): φ_R(r) = u_B(r)/u0(r), u_R(r) = u0(r) [exp(φ_R(r)) − 1].
5. V = {(x/b)² + (y/b)² + (z/a)² < 1}; ρ_e(r) := √((x/b)² + (y/b)² + (z/a)²) is the "elliptic radius" (1 on the surface S). Volume V = 4πab²/3.
6. Scalar ↔ vector for comparisons: u_gt := E_sca·conj(E0)/|E0|² with E0 the incident polarization vector of the sibling's `incident_plane_wave` (co-polarized projection). Cross-polarized energy is reported, never compared.

---

## 2. Mathematics

### 2.1 Reduction of the Born integral to a surface integral

Let w(r') := (k̂·r') u0(r') / (2ik). Then (∇'² + k²) w = u0, because
∇²[(k̂·r) u0] = 2ik u0 − k² (k̂·r) u0.
Green's second identity on V with outward normal n̂ (∫_V (G∇'²w − w∇'²G) dV' = ∮_S (G ∂_n w − w ∂_n G) dS') with ∇'²w = u0 − k²w and ∇'²G = −k²G − δ(r' − r) gives

∫_V G(|r−r'|) u0(r') d³r' = ∮_S [ G ∂_{n'}w − w ∂_{n'}G ] dS' − χ_V(r) w(r),

χ_V = 1 for r inside V, 0 outside (½ on S, not used — see 3.1). Hence

**u_B(r) = f { ∮_S [ G(R) ∂_{n'}w(r') − w(r') ∂_{n'}G(R) ] dS' − χ_V(r) w(r) }**, R = |r − r'|,

with ∇'w = k̂ u0(r') [1 + ik (k̂·r')]/(2ik) and ∂_{n'}G = (ik − 1/R) G(R) (n̂'·(r' − r))/R.

Remarks. (i) Exact for every r ∉ S, inside or outside; the interior 1/R singularity of the volume integrand never appears. (ii) u_B is a volume potential with bounded density, hence continuous across S; in the formula the double-layer term −∮ w ∂_nG jumps by ∓w across S and cancels the jump of χ_V w. An error in the jump term would show up as an O(1) error at interior points, which T1 checks against an independent series. (iii) The result does not depend on the choice of the particular solution w; T2 checks it against brute-force volume quadrature.

### 2.2 Surface quadrature

Parametrize S by r'(θ,φ) = (b sinθ cosφ, b sinθ sinφ, a cosθ), θ ∈ [0,π], φ ∈ [0,2π). Then
n̂' dS' = b sinθ (a sinθ cosφ, a sinθ sinφ, b cosθ) dθ dφ
(unnormalized normal times area element; every ∂_{n'}(·) dS' above is ∇'(·)·(n̂' dS')). The integrand is analytic and 2π-periodic in φ and analytic in θ on [0,π] (the sinθ factor removes the pole degeneracy), so
- θ: Gauss–Legendre with N_θ nodes on [0,π] (`gauss_legendre` mapped from [−1,1]);
- φ: N_φ = 2N_θ uniform nodes, weight 2π/N_φ (trapezoid, spectral for periodic analytic integrands);
converge spectrally for every r ∉ S. Nodes, weights, r', n̂'dS', u0(r'), w(r') and ∇'w(r') are precomputed once per N_θ and shared by all points.

**Kernel evaluation (cancellation-free).** Never form e^{ikR} from R directly: for |r| ≫ max(a,b) the rounding of R = |r − r'| (≈ eps·|r| per node) becomes a node-dependent phase error k·eps·|r| that the convergence check of 2.3 cannot see, because both quadratures share it (measured 2026-09-12: sphere, x = 100, |r| = 1e9 a, forward axis: N-check difference 2.7e-9 while the actual error was 2.0e-8, constant in N). Instead write, with ρ := |r|,
R − ρ = (|r'|² − 2 r·r')/(R + ρ)  (exact algebra, no cancellation; valid for every r including r = 0 since R + ρ > 0 on S),
G(R) = e^{ikρ} · e^{ik(R−ρ)}/(4πR), and likewise in ∂_{n'}G. The whole surface sum and the jump term (whose u0(r) = e^{ikρ}e^{ik(k̂·r − ρ)} is factored the same way) are accumulated *without* the common factor e^{ikρ}; the convergence check of 2.3 is applied to that phase-free quantity; e^{ikρ} is multiplied on as the very last operation. Consequence for tests: u_B·exp(−1i*k*ρ) with the same double product k*ρ recovers the phase-free field exactly (3.2). With this form the N-check difference at |r| = 1e9 a drops to 3e-15.

Start value: the phase of e^{ikR} u0(r') varies by at most ≈ 4 k max(a,b) over the surface, so N_θ0 = ⌈2 x⌉ + 16 with x = k max(a,b) (≈ 3 nodes per oscillation). Convergence slows for r close to S: the integrand then has a peak of width ≈ d = dist(r, S) and height ≈ 1/d², and N_θ ≳ 10 max(a,b)/d is needed (3.1).

### 2.3 Convergence loop (driven by `tol`)

For every point compute the phase-free Born values B_N and B_M with M = ⌈1.5 N_θ⌉ (N_φ = 2N_θ for each), then form the corresponding phase-free Rytov values R_N and R_M as in 2.4. Accept the finer values only if both
|B_M − B_N| ≤ tol · max(|B_M|, s(r)) and |R_M − R_N| ≤ tol · max(|R_M|, s(r)),
where s(r) := |f| V / (4π max(|r|, max(a,b))). This natural scale prevents pattern zeros from failing spuriously; `info.err` is the larger of the two normalized differences. Failed points are re-evaluated as a group with N_θ ← 2N_θ. **N_max rule:** ordinary doubling continues while its check grid fits. If the next doubled pair would exceed `N_max`, evaluate one final cap pair (⌊N_max/1.5⌋, N_max), so `N_max` is both the largest grid built and the finest available failure value. Example: `N_max = 200` with start 76 runs (76, 114), then (133, 200). An initial pair that does not fit is an input error. Remaining failures: `on_fail='error'` (default) raises an error naming the count and their min |ρ_e − 1| — nothing is returned in that case; `'warn'` returns the final values, `info.unconverged` true for them, `info.validated = false`. The same mask and policy apply to non-finite Born or Rytov output (3.3), so a usable finite field alongside failed points is available only in `'warn'` mode.

### 2.4 Rytov

φ_R = u_B ./ u0 (u0 has modulus 1 everywhere, no division hazard); u_R = u0 .* expm1(φ_R). `expm1` accepts complex input and avoids the cancellation of exp(φ)−1 for small φ (measured 2026-09-12 at φ = 1e-12(1+0.5i): `expm1` matches the series φ + φ²/2 to all digits, exp(φ)−1 is off by 8.9e-5 relative). Overflow: exp(φ_R) is infinite when Re φ_R > log(realmax) ≈ 709.8, which no physically meaningful input reaches; the output is then non-finite and handled by 3.3 — it is never silently reported as valid.

### 2.5 Reference solutions (tests only; live as local functions of the test script)

(a) **Sphere partial-wave series** (a = b), independent of 2.1–2.2. With the addition theorem G = ik Σ_l j_l(k r_<) h_l^{(1)}(k r_>) Σ_m Y_lm(r̂) Y*_lm(r̂') and e^{ik k̂·r'} = 4π Σ_l i^l j_l(kr') Σ_m Y*_lm(k̂) Y_lm(r̂'), orthonormality and Σ_m Y Y* = (2l+1)/(4π) P_l give
u_B(r) = f · ik · Σ_{l≥0} (2l+1) i^l P_l(k̂·r̂) I_l(r),
I_l(r) = ∫_0^a r'² j_l(kr') j_l(k min(r,r')) h_l^{(1)}(k max(r,r')) dr',
i.e. for r ≥ a: h_l(kr) ∫_0^a r'² j_l² dr'; for r < a: h_l(kr) ∫_0^r r'² j_l² dr' + j_l(kr) ∫_r^a r'² j_l h_l dr'.
Radial integrals with `integral` (RelTol 1e-12, AbsTol 0; smooth integrands). j_l, y_l = √(π/2z) J_{l+½}, Y_{l+½}; P_l by the three-term recurrence. l_max = ⌈ka + 4.05 (ka)^{1/3} + 20⌉, verified by adding 10 terms (change < 1e-12). r = 0: only l = 0 (j_0(0) = 1). Sanity: for ka → 0, u_B(0) → f a²/2.
Far points (r ≥ a, kr ≫ l_max²): evaluate h_l with the common phase removed using the exact finite sum (DLMF 10.49.6)
h_l^{(1)}(z) = e^{iz} Σ_{m=0}^{l} i^{m−l−1} (l+m)! / (m! (l−m)! 2^m z^{m+1}),
i.e. h̃_l(z) := h_l e^{−iz} is a polynomial in 1/z (well conditioned for z ≫ l²; used only when z > 10 l_max², otherwise Bessel routines). Then u_B(r) e^{−ikr} = f ik Σ_l (2l+1) i^l P_l(k̂·r̂) h̃_l(kr) ∫_0^a r'² j_l² dr' is an independent, phase-free far-field reference that does not rely on the far-field approximation, matching the factored output of 2.2.
(b) **Brute-force volume quadrature** (any a, b; exterior points only): r' = (b s sinθ' cosφ', b s sinθ' sinφ', a s cosθ'), dV = ab² s² sinθ' ds dθ' dφ'; Gauss–Legendre in s ∈ [0,1] and θ' ∈ [0,π], uniform φ'. Analytic integrand for r ∉ V̄ ⇒ spectral; use N = 144 per dimension for x ≤ 10 after verifying N = 96 versus 144 within 1e-10. The originally planned 64/96 pair missed that threshold for the prolate case.
(c) **Rayleigh–Gans far field** (any a, b): with |r − r'| ≈ r − r̂·r' in G,
u_B(r) → f V e^{ikr}/(4πr) · F(q), q = k (r̂ − k̂), F = 3 (sin Q − Q cos Q)/Q³, Q = √(b² (q_x² + q_y²) + a² q_z²), with F → 1 − Q²/10 for Q < 1e-3.

---

## 3. Numerical safeguards

### 3.1 Points near or on the surface
`inside = ρ_e < 1` decides χ_V. Points with |ρ_e − 1| < 1e-9 are rejected with an error ("offset the point from the surface"): on S the double-layer integrand is not integrable by the plain rule and the ½-jump formula would need a principal value. Points with small but nonzero distance are handled by the convergence loop; the tested envelope (Section 7) uses the conservative physical-clearance proxy |ρ_e − 1| min(a,b) ≥ 0.02 max(a,b). The unscaled elliptic-radius condition is insufficient at aspect-ratio extremes. `ponytail:` global-N refinement per failing group; add local θ–φ subdivision around the foot point of r on S if near-surface fields are ever needed.

### 3.2 Far points
The kernel of 2.2 keeps every node-dependent phase cancellation-free, so the phase-free field is subject to the error criterion of 2.3 and to the independent far-point references of 2.5(a) (T1 far cases) and 2.5(c) (T3). What remains is the common factor e^{ikρ}: in this double implementation the product k·ρ carries a rounding error ≈ k ρ eps, e.g. ≈ 2e-5 rad at ρ = 1e9 max(a,b), x = 100. That is a property of this implementation, not a fundamental limit: treating the double inputs as exact, k ρ mod 2π could be reduced in higher precision; this is not done. The returned u0 is reconstructed with the same factored phase as u_B and u_R, so the advertised Rytov identity remains internally consistent even at far off-axis points. Accuracy statement, therefore: |u_B|, |u_R| and the fields with the common phase removed obey `tol`; the absolute phase of u0, u_B and u_R carries an additional error ≈ k ρ eps. Tested envelope: |r| ≤ 1e9 max(a,b) (T1 far cases go to 1e9 a). No far-field kernel is implemented (`ponytail:` replace G by e^{ikρ}e^{−ik r̂·r'}/(4πρ) if |r| ≫ 1e9 max(a,b) is ever needed). Note that the RGD reference of 2.5(c) cannot verify at the `tol` level: its 1/r truncation was 2.0e-8 for the x = 100 sphere on the forward axis at r = 1e9 a (measured 2026-09-12), and is O(k max(a,b)²/r) in general; T3 therefore uses r = 1e7 max(a,b), x ≤ 8 and tolerance 1e-5, while tol-level far checks are T1's job.

### 3.3 Overflow / underflow
`expm1(φ_R)` overflows when Re φ_R exceeds ≈ 709.8 (2.4), and finite inputs can also overflow intermediate or Born arithmetic at extreme magnitudes. Both quadrature pairs and final `uB(:)` and `uR(:)` are therefore checked with `isfinite`; any non-finite point is added to `info.unconverged`, makes `info.validated` false, and follows `opts.on_fail` exactly as in 2.3 (error by default; `'warn'` returns both arrays with the mask).

### 3.4 Cost and memory
Per point 2 N_θ² kernel evaluations (× 2.25 for the 1.5 N_θ check). x = 30 ⇒ N_θ = 76, ≈ 2.6e4 nodes ⇒ 1e4 points ≈ 3e8 complex exponentials. Points are processed in chunks so that (points per chunk) × (nodes) ≤ `opts.chunk`. The default is 2e6 point-node entries: one complex matrix is ≈ 32 MB, while the several simultaneous real and complex kernel matrices keep the practical working set in the few-hundred-MB range. The earlier 2e7 default could require 1–2 GB and is not used. No `parfor` (`ponytail:` chunks are the place to add it).

---

## 4. Files and interfaces (`forward/`)

| File | Purpose | Interface |
|---|---|---|
| `born_rytov_spheroid.m` | Entry point: validate inputs and finite derived k/f, k̂; surface nodes (2.2); chunked evaluation of 2.1; convergence loop (2.3); Rytov (2.4); reshape. Local functions `surface_nodes(a,b,k,khat,N_theta)` and `born_eval(nodes,pts,f,k,khat,inside)`. | `[uB, uR, u0, info] = born_rytov_spheroid(lambda, n_p, n_m, a, b, theta_inc, phi_inc, X, Y, Z, opts)`. `opts` fields: `tol` (1e-8), `N_theta` (start value; default ⌈2x⌉ + 16 from 2.2), `N_max` (1200, largest grid built — rule in 2.3), `on_fail` ('error' \| 'warn'), `chunk` (2e6 point-node entries). Outputs have the size of `X`. `info`: `k`, `f`, `khat`, `V`, `N_theta` (per point, after refinement), `err` (maximum normalized Born/Rytov check), `inside`, `unconverged` (convergence failures and non-finite output), `validated`, `time`. |
| `compare_with_spheroid.m` | Section 6: sweeps, observables, tables and figures against `mie_field` / `spheroid_field`. | script; prints the tables that are copied into Section 9. |
| `tests/test_born_rytov.m` | T1–T5; local functions `born_sphere_series`, `born_volume_quadrature`, `rgd_farfield` (2.5). | script with asserts |
| `tests/run_all.m` | Adds the sibling folder to the path, runs `test_born_rytov`, prints pass/fail and timing. | script |

Everything user-facing is in the user's length unit (shared by `lambda, a, b, X, Y, Z`), as in the sibling.

---

## 5. Tests (definition of done) — `tests/test_born_rytov.m`

T1 **Sphere series** (checks 2.1 including the jump term, 2.2, and the convergence loop). a = b; cases (x, n_p/n_m) ∈ {(0.5, 1.05), (5, 1.01), (30, 1.02+0.01i)}; ζ = 40°, φ_i = 20°. Points: 50 random interior (ρ_e ≤ 0.9), 50 random exterior (1.1 ≤ ρ_e ≤ 10), the origin, 5 points on the z-axis inside and outside. Require |u_B − u_series| ≤ 1e-8 max|u_series| and `info.validated`.
Far cases (checks the cancellation-free kernel of 2.2 against the phase-free series of 2.5(a)): x ∈ {5, 100}, n_p/n_m = 1.01, ζ = 0; 20 random directions at each |r| ∈ {1e4, 1e7, 1e9} a, including the forward axis. Compare the phase-free quantities u_B·exp(−1i*k*|r|) (same double product as inside the solver) with the h̃_l series; require agreement ≤ 1e-8 relative. Also assert that the *naive* kernel e^{ikR} (kept as a test-only local function) fails this at |r| = 1e9 a, x = 100 (expected ≈ 2e-8), so the test guards the reason for 2.2.

T2 **Volume quadrature** (checks the surface reduction for a ≠ b). (a/b, x) ∈ {(2, 8), (0.5, 8)}, n_p/n_m = 1.1, ζ = 30°, φ_i = 0; 30 random exterior points with 1.2 ≤ ρ_e ≤ 5. Require |u_B − u_vol| ≤ 1e-8 max|u_vol| after the volume rule itself has been verified (N = 96 versus 144 within 1e-10); use the N = 144 value as the reference.

T3 **Rayleigh–Gans far field** (checks geometry a ≠ b and incidence handling at large |r|). Same two spheroids; r = 1e7 max(a,b); 181 directions in the plane spanned by k̂ and ẑ (θ from 0 to 180° relative to k̂) plus 50 random directions. Require |u_B r e^{−ikr} − f V F(q)/(4π)| ≤ 1e-5 |f| V/(4π).

T4 **Rytov**. (i) n_p = n_m ⇒ u_B = u_R = 0 exactly (f = 0, no NaN). (ii) u0 + u_R = u0 exp(u_B/u0) to 1e-14 relative at points where |φ_R| > 1e-3. (iii) Small-phase accuracy: for x = 1, n_p/n_m = 1 + 1e-9 (|φ_R| ≈ 1e-9 … 1e-12) at 20 exterior points, |u_R/u0 − (φ_R + φ_R²/2)| ≤ 1e-8 |φ_R| (fails with exp(φ)−1, passes with `expm1`). (iv) Born and Rytov agree to first order (Kak p. 217): for x = 5, n_p/n_m = 1 + 1e-4, ζ = 0, at 100 random points with 1.2 ≤ ρ_e ≤ 3: |u_R − u_B| ≤ 1e-3 |u_B|. (v) Non-finite output: with a deliberately absurd f (n_p/n_m = 1 + 1e6, x = 30, an exterior point) u_R is not finite; default `on_fail` errors, `'warn'` returns finite u_B, `info.unconverged` true at that point and `validated` false. (vi) Rytov-sensitive convergence: use a finite case whose initial pair passes the Born criterion but fails the transformed Rytov criterion; require refinement and agreement with an independently converged sphere-series Rytov reference within `tol`. (vii) At an off-axis point with |r| = 1e9 a and finite nonlinear φ_R, require u_R = u0 expm1(u_B/u0) to 1e-12 relative using the returned fields; this catches inconsistent common-phase evaluation.

T5 **Convergence, envelope, failure modes**. (i) Over the corners of Section 7 (x ∈ {1, 30, 100}, a/b ∈ {1/4, 1, 4}, ζ ∈ {0, 60°}), 20 random points each satisfying |ρ_e − 1| min(a,b) ≥ 0.02 max(a,b): `info.validated` true, and re-running with `opts.N_theta` doubled changes u_B by ≤ 1e-8 relative. (ii) A sphere point at ρ_e = 1.02 converges within the default `N_max`. (iii-a) Surface rejection: a point at ρ_e = 1 + 1e-10 raises the "offset the point from the surface" error regardless of `on_fail`. (iii-b) Convergence failure, fully specified: sphere a = b, x = 30, n_p/n_m = 1.1, ζ = 0, the single point (1 + 1e-4) a on the +x axis, `opts.N_theta = 76`, `opts.N_max = 200`. By the N_max rule of 2.3 the loop evaluates (76, 114), then the final cap pair (133, 200), and remains above `tol`. Assert: default `on_fail` raises the convergence error; `'warn'` returns `info.unconverged` true, `validated` false, `info.N_theta = 200`, the N = 200 value, and `info.err` equal to the maximum normalized Born/Rytov difference of the final pair. A neighboring-limit regression with `N_max = 115` must evaluate (76, 115), not stop at 114. If either point converges instead, the test fails explicitly. (iv) Reject an initial check pair above `N_max` and finite inputs that produce non-finite derived k/f. (v) Reshape: X of size 7×3×2 gives outputs of that size; scalar inputs work.

Coverage of the ground-truth comparison is Section 6 (mostly a report, one pass/fail item).

---

## 6. Validation against the vector ground truth — `compare_with_spheroid.m`

Ground truth: every case is routed through `spheroid_field`, including a = b (for which the sibling routes internally to `mie_field`), with `field='scattered'`, `tol=1e-6`, `on_fail='error'`, and `pol=[1,0]` (E0 = ê1 of sibling §1.3, |E0| = 1); u_gt = E_sca·conj(E0) (1.6). Cross-polarized fraction ‖E_sca − u_gt E0‖/‖E_sca‖ is printed alongside. Sections 1–5 do not depend on these functions.
Status (2026-09-12): both sibling APIs exist with the planned signatures (`[E,H] = mie_field(lambda,n_p,n_m,a,theta_inc,phi_inc,...)`, `[Ex,Ey,Ez,info] = spheroid_field(lambda,n_p,n_m,a,b,...)`). "Exists" is not "validated as ground truth": comparison results are accepted only after the sibling's complete `tests/run_all.m` is green. Within the comparison, non-finite fields or `info.validated=false` are eligibility errors, not additional model-accuracy assertions; the sole accuracy assertion remains the item below. Compact validation metadata is retained, but the unused full modal solution in `info.sol` is discarded before cases are accumulated.

What the scalar model can and cannot match (state this in the script header):
- The scalar Born amplitude differs from the *vector first-order far-field* amplitude by the dipole factor (I − r̂r̂)·ê0, i.e. by 1 − (r̂·ê0)² in the co-polarized projection. In the H-plane (plane containing k̂ and ê2, i.e. r̂ ⊥ ê0) this factor is exactly 1. In the E-plane it is not, and the vector first-order term vanishes at 90° while the exact finite-contrast vector field need not.
- In the near and interior field the scalar/vector difference is O(1) in n_δ but small for k·(object size) ≫ 1 (Kak p. 209: depolarization negligible when the wavelength is much smaller than the inhomogeneity). Small-x interior comparisons therefore measure the scalar approximation, not Born.

Observables (for each case: Born, Rytov; relative L2 error ‖u_X − u_gt‖/‖u_gt‖):
- O1 Detector plane ⊥ k̂ through r_D = 2 max(a,b) k̂, spanned by ê1, ê2, half-width 3 max(a,b), 101×101 points. Errors on the whole plane, and separately along the ê2-line (H) and ê1-line (E). Rytov-domain error ‖φ_X − φ_gt‖/‖φ_gt‖ with φ_gt = log(1 + u_gt/u0) is computed on the two lines only (H-line, E-line), where 1-D unwrapping is well defined: take the principal log, then `unwrap` the imaginary part starting from the outer end of the line (where the scattered field is weakest) inward. The principal-value check |Im φ| < π proves nothing (the principal log always satisfies it), so a line is accepted for this metric only if all three practical checks pass: (a) at the reference end |u_gt/u0| < 1e-2, so φ ≈ 0 there fixes the branch; (b) along the whole line |1 + u_gt/u0| > 0.1, so the log is never near its singularity; (c) recomputing the ground truth with half the sample spacing and unwrapping again reproduces the phase at the original samples to 1e-3 rad. These are practical criteria that catch the aliasing seen in practice, not a proof that no aliasing is possible. Lines failing any check are marked "skip" in the table with the failing check named.
- O2 Far-field rings, r = 1e7 max(a,b), θ ∈ [0°, 180°] step 1°: H-plane and E-plane. In the E-plane additionally show u_X·(1 − (r̂·ê0)²) to demonstrate that the residual there is the scalar/vector effect.
- O3 Axial line along k̂ through the centre, from −3 to +3 max(a,b), including interior points that satisfy |ρ_e − 1| min(a,b) ≥ 0.02 max(a,b): Re/Im of u0 + u_X versus u0 + u_gt and the accumulated phase.

Sweep: n_δ = (n_p − n_m)/n_m ∈ {1e-3, 1e-2, 3e-2, 1e-1} × x ∈ {1, 5, 20}; shapes sphere, prolate a/b = 2, oblate a/b = 1/2; ζ ∈ {0, 30°} (φ_i = 0). Report the accumulated phase Δφ = 2 k L_k n_δ (L_k = spheroid half-length along k̂) next to every row: Born is expected to hold while Δφ ≲ π/2 (Kak 50–51: a n_δ < λ/4), Rytov further but limited by the index gradient at the boundary (Müller 3.47–3.52), and both degrade at small x for the scalar reasons above.

Pass/fail item (the only assertion): sphere, x = 5, ζ = 0, O2 H-plane, θ ∈ [0°, 60°]: at n_δ = 1e-3 the relative error is ≤ 2e-2 for both u_B and u_R (expected ≈ x n_δ = 5e-3), and error(n_δ = 2e-3)/error(n_δ = 1e-3) ∈ [1.5, 2.5] (first-order scaling of the model error).

---

## 7. Tested envelope

| Quantity | Range covered by T1–T5 |
|---|---|
| x = k_m max(a,b) | 0.5–100 (T1 covers both endpoints; T2–T4 use 1–30) |
| a/b | 1/4 … 4 |
| n_p/n_m | Finite complex input is accepted only when derived f is finite and output passes the failure policy; T1 tests 1.02+0.01i, while Section 6 tests real 1.001–1.1. Extreme finite values are numerical failure-path cases in T4, not part of the validated physical envelope. |
| Evaluation points | |ρ_e − 1| min(a,b) ≥ 0.02 max(a,b), |r| ≤ 1e9 max(a,b) (phase-free field to `tol`; absolute phase carries ≈ k|r|·eps, 3.2); inside and outside |
| ζ | any; φ_i any |
Outside: solver still runs; near-surface points may exhaust `N_max` (3.1).

---

## 8. Notes and deliberate simplifications

- **Vector Rytov** (asked during planning): this implementation is restricted to the scalar model by scope — the two references are scalar, and the ground-truth comparison uses the co-polarized projection (1.6). Vector Rytov formulations exist: a component-wise log E_j = E0_j exp(φ_j) is the naive one (undefined where E0_j = 0), and a well-posed version applies the Rytov ansatz to the scattering (Jones) matrix instead of to field components: Oh et al., "Extending Rytov approximation to vector waves for tomography of anisotropic materials", Phys. Rev. Lett. 134, 068101 (2025), arXiv:2404.17206. That, or a vector Born with the dyadic Green's function Ḡ = (I + ∇∇/k²)G, is the extension path; the surface reduction of 2.1 applies to the vector Born volume integral term by term.
- **Direct summation in practice** (asked during planning): for voxel objects the production route is FFT convolution on the grid (Lippmann–Schwinger / DDA solvers) or the Fourier diffraction theorem from a plane (Müller §5.1), both O(N log N); direct O(N_vox · N_pts) sums are used for validation or few arbitrary points. For a homogeneous spheroid the exact surface reduction of 2.1 is cheaper (2-D instead of 3-D nodes) and more accurate (no staircase, spectral) than either, which is why this plan uses it.
- Deliberate simplifications, each marked `ponytail:` in code: no local refinement near S (3.1); no far-field kernel and no high-precision reduction of the common phase (3.2); no `parfor` (3.4). Plane-wave illumination only: a focused beam would need a particular solution w of (∇² + k²)w = u0 — available for a plane-wave spectrum, otherwise fall back to 2.5(b)-style volume quadrature.

---

## 9. Implementation order (each step ends with its verification)

1. `born_rytov_spheroid`: nodes (2.2), kernel and jump term (2.1), fixed `N_theta`, no loop → T2 (volume quadrature) and T3 (RGD far field), both against reference functions written first.
2. Convergence loop (2.3), chunking (3.4), `opts`, `info`, surface rejection (3.1) → T1 (sphere series, interior points included) and T5.
3. Rytov line (2.4) → T4. `tests/run_all.m`.
4. `compare_with_spheroid.m` through `spheroid_field` for all shapes (with its internal sphere-to-Mie route) after the sibling's complete suite is green → the pass/fail item, 72-case table, and nine figures of Section 6. Append the obtained summaries as a results note at the end of this file.
