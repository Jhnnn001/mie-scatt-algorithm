# Where first Rytov fails for this prolate

Date: 2026-09-12. All new work is confined to `forward/exp`.
Neither the production forward model nor the sibling electromagnetic solver was changed.
The legacy `forward/born_rytov_spheroid.m` was not used.

## Conclusion

The residual error is predominantly first-Rytov truncation, not scalar reduction or voxel resolution.
The relevant neglected term is the square of the **complex logarithmic field gradient**;
diffraction itself is retained by first Rytov.
There is a numerically verified zero of the exact scalar total field immediately next to
the chosen detector plane, where that logarithm cannot be smooth or finite.

This is not evidence that homogeneous spheroids are intrinsically harder than real cells.
It demonstrates that this geometry, contrast, wavelength, and observation plane are not
a uniformly weak-log-gradient configuration, despite their simple material description.

## Configuration and error metric

| Quantity | Value |
|---|---:|
| Vacuum wavelength | 0.532 um |
| Medium index | 1.335381534 |
| Particle index | 1.37, real |
| Prolate semiaxes | a = 5 um along z; b = 2.5 um transversely |
| Incidence / electromagnetic reference polarization | +z / x |
| Detector NA | 0.1 |
| Detector planes | z = 5.2, 7, 14 um; object tip at z = 5 um |
| Estimated central ray phase delay | 4.088613 rad |
| Continuous-shape spectral window | 51.2 um; 512 or 1024 lateral samples |

Reported detector errors are `norm(U_s_approx - U_s_exact) / norm(U_s_exact)`
for the **complex scattered field**, after the same NA pupil.
They are not intensity errors, relative total-field errors, or reconstruction errors.
The shape transform is analytic here: no voxelized object or Ewald interpolation is needed.

## 1. Separate scalar reduction from Rytov truncation

An experimental scalar reference solves the scalar Helmholtz boundary problem using
spheroidal angular and radial functions already present in `spheroid-analytic-forward`.
It enforces continuity of U and its normal derivative, not Maxwell boundary conditions.
It is restricted to real-index prolates with axial scalar plane-wave incidence.

At NA = 0.1, exact scalar scattering differs from the vector Ex reference by **0.252653%**.
First Rytov at z = 7 um differs from exact scalar by **24.1742%**, and from vector Ex by
24.3021%, so scalar reduction is not the principal explanation for the earlier error.

Independent/reference checks:

| Check | Result |
|---|---:|
| Zero contrast: total plane wave and gradient | PASS |
| Near-sphere comparison against independent scalar spherical partial waves | 1.12e-6 relative field difference; shape mismatch 1e-6 |
| Physical prolate boundary residual at independent collocation nodes | 5.47e-15 |
| Scalar far amplitude, L = 112 versus 132 | 1.21e-12 relative difference |
| Analytic far limit versus finite-radius field | 7.51e-6 at r = 5e6 um |
| Exact weak-contrast derivative versus analytic Born spectrum | 5.10e-5 |
| Second perturbation coefficient, contrast step 0.01 versus 0.005 | 8.36e-5 |
| N = 512 versus 1024, largest change in first-Rytov error metric | 1.28e-9 |

Order agreement is a convergence check, not a universal absolute-accuracy guarantee.
The scalar-vs-vector comparison uses the sibling validated vector solver and finite-radius
far extraction; the scalar far amplitude uses its analytic outgoing asymptotic limit.

## 2. The term actually omitted

For U0 = exp(ikz), define ψ = log(U/U0) and f = k0²(n² - nm²).
Where U is nonzero, the exact scalar equation is

`∇²ψ + 2ik ∂zψ + Q = -f`

`Q = (∂xψ)² + (∂yψ)² + (∂zψ)²`.

There is **no complex conjugation** in Q.
First Rytov drops Q but keeps ∇²ψ, which is why it includes diffraction.
The nonlinearity is in the logarithmic representation, not nonlinear optical material response.
This distinction is consistent with the diffraction-aware Rytov treatment in
[Sung et al., 2009](https://web.mit.edu/spectroscopy/doc/papers/2009/optical_sung_09.pdf).

We evaluate `∂jψ = (∂jU)/U - i kj` from the scalar reference and analytic field derivatives:
no principal-log branch or phase-unwrapping assumption enters this diagnostic.
Independent finite differences at interior/exterior points verify the exact logarithmic
equation to a normalized residual of 4.57e-7 for the physical contrast.

| Interior region | Volume-weighted RMS magnitude of Q / f, grid 0.05 um | Grid 0.025 um |
|---|---:|---:|
| Central cylinder, rho < b/2, restricted to the object | 0.2410 | 0.2394 |
| Outer region, rho > 0.8 b, restricted to the object | 0.2863 | 0.2858 |

These are regional RMS ratios, not a pointwise bound or a formula for detector error.
The central region is heterogeneous in this diagnostic: its unweighted median |Q|/f
is only about 0.019, whereas the outer-region median is about 0.200.
The transverse contribution accounts for norm ratios 0.915 and 0.990 respectively;
these are complex-term norm ratios, not energy fractions.
At one tenth of the scattering potential, the RMS ratios fall to 0.0149 and 0.0330.

## 3. A direct failure point near z = 7 um

Starting from the fine spatial scan, a two-variable Newton solve finds a zero of U at

`rho = 0.1598748398 um, z = 6.996322478 um`.

The residual |U| is 2.41e-14 at L = 112 and 8.64e-12 when independently reevaluated at
L = 132; the incident amplitude is one.
Axisymmetry makes this an azimuthal ring of zeros in the scalar reference.
Consequently ψ is singular there and its gradient cannot satisfy a uniformly small-gradient
assumption in a neighborhood of that ring.
The zero is about 0.00368 um before the z = 7 um detector plane.

This is the **pre-objective** field, not an assertion that the NA-filtered measured field
has the same zero, nor proof that this one ring alone accounts for all 24.2% detector error.
The interior Q calculation and the contrast controls provide separate evidence.
Exterior f is zero, so a ratio Q/f is not used there.
Exterior RMS statistics in the CSV exclude |U| <= 0.1 and remain sensitive to spatial
sampling near zeros; they must not be treated as converged global validity bounds.

## 4. Controlled changes

### Observation distance and perturbation order

| z (um) | First Rytov error | Second Rytov diagnostic error |
|---|---:|---:|
| 5.2 | 16.8730% | 4.3934% |
| 7 | 24.1742% | 13.3812% |
| 14 | 47.2136% | 2704.996% |

The second coefficient is extracted from independent scalar solutions at small positive
and negative scattering-potential scales, not fitted to the physical-contrast target.
With U(t) = U0 + t U1 + t² U2 + ..., it is
`ψ1 = U1/U0`, `ψ2 = U2/U0 - ψ1²/2`; the experiment evaluates `U0 exp(t ψ1 + t² ψ2)`.
The logarithmic equation gives `L ψ2 = -∇ψ1·∇ψ1`, with `L = ∇² + 2ik ∂z`.
The substantial second-order effect supports the truncation diagnosis, but the z = 14
failure explicitly rules out recommending second order as a generally reliable fix.

For free propagation, the exact total field obeys the linear Helmholtz equation; its
logarithm still obeys the equation containing Q even in uniform water.
Thus applying first Rytov separately at successive detector planes need not agree with
forming an approximate field once and propagating that field exactly.
This is not the assertion that diffraction is missing from first Rytov.

### Contrast, at z = 7 um

The parameter t multiplies f, with `n(t) = sqrt(nm² + t (np² - nm²))`.

| t | First Rytov error | Second Rytov error |
|---|---:|---:|
| 0.10 | 2.5164% | 0.0459% |
| 0.25 | 6.1535% | 0.2651% |
| 0.50 | 12.0054% | 1.2326% |
| 1.00 | 24.1742% | 13.3812% |

### Geometry, preserving the central phase delay

The wavelength, index contrast, 10 um axial thickness, detector z = 7 um, and NA are fixed.

| Geometry | First-Rytov scattered-field error |
|---|---:|
| Original prolate: transverse radius 2.5 um | 24.1742% |
| Sphere: transverse radius 5 um | 12.7176% |
| Infinite uniform slab: thickness 10 um | 2.9774% |

The slab reference includes both interfaces and internal reflections; its relative
total-field error is 5.3004%, a different denominator from the scattered-field table.
The sphere uses an independent exact scalar spherical partial-wave series, not a
nearly spherical focal-coordinate limit.
These controls show that a central phase delay of 4.09 rad alone does not determine
Rytov accuracy, and that a broader transverse geometry improves this particular case.

The ray estimate `phase(rho) ≈ k0 (np - nm) thickness(rho)` offers intuition:
the same phase difference spread across a shorter lateral distance produces a steeper
wavefront, even if the thickness is continuous.
For a spheroid `thickness(rho) = 2a sqrt(1 - rho²/b²)`; continuity does not bound its slope,
and it says still less about the exact diffracted field near destructive-interference zeros.
This ray argument is interpretation, not the reference solution used in the experiments.

## Scope and reproducibility

The diagnostic is complete for the specified axial prolate and narrow-NA field metric;
it does not establish comparative accuracy for actual biological cells or inverse reconstructions.
The spectral experiments retain the production model's forward-propagating hemisphere
and finite transverse window; they omit evanescent modes.
The direct spatial Q and field-zero calculations do not use that spectral truncation.
The prior `forward/tests/rytov_analysis.csv` additionally documents window-size and
voxel controls; its 128-grid production error is 25.735%, versus 24.174% in this
continuous-shape scalar-reference diagnostic, with different numerical sampling.

Run from the repository root in MATLAB:

```matlab
addpath('forward/exp', 'spheroid-analytic-forward');
test_scalar_prolate;
diagnose_rytov(512);
diagnose_rytov(1024);
diagnose_rytov_term(.05);
diagnose_rytov_term(.025);
diagnose_rytov_geometry;
locate_field_zero;
addpath('forward/tests');
test_born_rytov_voxel;
```

`verification.log` records the final rerun; `diagnosis.log` preserves preliminary probes,
including the rejected near-sphere ODE route.
CSV files contain comparison/validation metrics and the field-zero location; MAT files
contain spectral arrays and spatial maps for later plotting.
Redundant experiment-generated full ODE caches were removed; scripts regenerate the solver.
No new dependency, production second-order model, commit, or push was added.
