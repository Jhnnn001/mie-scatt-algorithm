# Oblate spheroid: xz total-Ex maps and Rytov error diagnosis

Computed September 12, 2026, using MATLAB R2024a Update 9; plotted with Matplotlib.

The requested figure is [oblate_Ex_xz_4x3.png](oblate_Ex_xz_4x3.png), with a
[PDF copy](oblate_Ex_xz_4x3.pdf). Rows are incidence angles 0, 5, 10, and 15 degrees;
columns are Exact, Born, and first Rytov. All panels share the color range 0–2.8.
The main numerical conclusion is that second-order Rytov reduces the scalar
reference detector error from 9.89–10.48% to 3.85–4.14% for these four angles.
It explains a substantial part of the first-order discrepancy, without eliminating it.

## What is plotted

- Homogeneous oblate spheroid: semiaxes b=5 um along x/y and a=3 um along z;
  full dimensions 10 x 10 x 6 um, background index 1.335381534, object index 1.365.
- Vacuum wavelength 0.532 um; incidence in the xz plane, phi=0, toward positive z.
- Original `pol=[1,0]` convention: k-hat=(sin(theta),0,cos(theta)) and
  E0=(cos(theta),0,-sin(theta)); the incident vector has unit magnitude.
- Color is total |Ex|/|E0|, including the incident field, on y=0. It is neither
  intensity |Ex|^2 nor scattered-field magnitude. The outline is the exact cross section.
- Raw fields are evaluated inside and outside the object, without a detector pupil.
  The grid is x=-8 to 8 um, z=-5 to 9 um, with 0.08-um spacing in both directions.
  The displayed arrows give the incident propagation direction.

Exact means the converged Maxwell spheroidal expansion. Born and Rytov use the
scalar Helmholtz model for the same continuous homogeneous spheroid, projected
onto laboratory x by cos(theta). These panels evaluate the full free-space Green
function; they are not images from the older voxel geometry or reconstructions
of refractive index. The existing voxel function evaluates a detector plane,
so its forward-only spectral construction is not used inside the object.

Let U0 be the scalar incident plane wave, and UB the first Born *scattered* field.
The plotted model fields are

$$E_{x,B}=\cos\theta\,(U_0+U_B).$$

$$E_{x,R}=\cos\theta\,U_0\exp(U_B/U_0).$$

## Errors associated with this figure

The amplitude error is norm(|Ex_model|-|Ex_exact|)/norm(|Ex_exact|), over the
entire displayed xz window. Complex error instead uses norm(Ex_model-Ex_exact)
in the numerator and therefore includes phase disagreement. All table entries
are percentages; these are total-field errors and depend on the stated window.

| Angle (deg) | Born amplitude | Rytov amplitude | Rytov complex | Scalar exact vs Maxwell, complex |
|---:|---:|---:|---:|---:|
| 0 | 58.689 | 5.220 | 8.275 | 2.242 |
| 5 | 58.928 | 5.291 | 8.361 | 2.276 |
| 10 | 59.695 | 5.494 | 8.601 | 2.437 |
| 15 | 60.962 | 5.808 | 9.000 | 2.687 |

Born substantially overestimates the downstream amplitude. At normal incidence,
the maximum |Ex| in the window is 1.470 for Exact, 2.774 for Born, and 1.590 for
Rytov. A centerline ray accumulates a phase delay of 2.099 rad through the object;
this is not a small phase for the linear total-field Born approximation.
Rytov retains an exponential phase response, but approximates its logarithm.
The raw scalar-versus-vector error above is distinct from the earlier 0.70%
co-polarized scattered-field error after a detector pupil.

## Why first Rytov still has about 10% detector error

The earlier 10.509% value was an illumination-pupil-area-weighted mean over
81 directions up to illumination NA=0.5, against a co-polarized Maxwell reference.
It used the complex scattered field at z=5 um after detector NA=1.0.
The four angles in this figure are a subset of that angular range, not a replacement
for the 81-direction average. Their maximum illumination NA is approximately 0.346.

To isolate the scalar forward approximation, the following new diagnostic uses
an exact scalar boundary-value reference at the same detector plane and NA.
Its error is norm(P[Us_model-Us_exact])/norm(P[Us_exact]), where P is the ideal
NA=1 circular pupil. It is not a percent error in the colors shown above.

| Angle (deg) | Born | First Rytov | Second Rytov diagnostic |
|---:|---:|---:|---:|
| 0 | 89.047 | 9.890 | 3.847 |
| 5 | 89.273 | 9.957 | 3.878 |
| 10 | 89.968 | 10.149 | 3.960 |
| 15 | 91.139 | 10.479 | 4.138 |

For scalar Helmholtz scattering, define f=k0^2(n^2-nm^2), km=k0 nm, and
U=U0 exp(psi), with U0=exp(i km k-hat dot r). Substitution gives

$$\nabla^2\psi+2i k_m\hat{\mathbf k}\cdot\nabla\psi+Q=-f,$$

$$Q=\nabla\psi\cdot\nabla\psi.$$

The dot product in Q has **no complex conjugation**. First Rytov omits Q and
sets psi1=UB/U0; it still retains the Laplacian and therefore diffraction.
Small index contrast alone does not make this omitted term vanish, and validity
of the scalar approximation does not establish validity of first Rytov.
See [Müller, Schürmann, and Guck, The Theory of Diffraction Tomography,
section 3.2, equations 3.34–3.38](https://arxiv.org/pdf/1507.00466).

The exact scalar field gives an interior xz-slice median |Q|/f of 0.038–0.039,
90th percentile 0.194–0.227, and RMS 0.466–0.496 across the four angles.
The term is spatially nonuniform; the 99th percentile is 1.69–1.77.
No near-zero field values required masking in this slice (minimum interior
|U| is 0.547 or larger). These are local two-dimensional diagnostics, not volume
averages or a formula for predicting the detector error percentage.

For a quantitative control, scale the scattering potential by t while preserving
shape and background, using n(t)=sqrt(nm^2+t(np^2-nm^2)). Expand the exact scalar
field as U(t)=U0+t U1+t^2 U2+..., then define

$$\psi_1=U_1/U_0,$$

$$\psi_2=U_2/U_0-\tfrac12(U_1/U_0)^2.$$

The second diagnostic uses U0 exp(t psi1+t^2 psi2). Its coefficients are extracted
from symmetric weak-contrast solves, not fitted to the physical-contrast result.
For consistency with the detector forward algorithm, this diagnostic forms psi1
and psi2 from forward propagating spectra before applying the pupil; the xz maps
and local Q above instead use the full Green-function fields. Second Rytov is
an error-diagnosis experiment, not a modification of the production forward function.

At normal incidence, the controlled contrast sweep gives:

| Potential scale t | Object index | First Rytov error (%) | Second Rytov error (%) |
|---:|---:|---:|---:|
| 0.1 | 1.338372877 | 0.867 | 0.034 |
| 0.25 | 1.342847397 | 2.193 | 0.214 |
| 0.5 | 1.350271980 | 4.512 | 0.880 |
| 1 | 1.365000000 | 9.890 | 3.847 |

The decrease with contrast and the improvement from the independently extracted
second term support truncation of the logarithmic-field expansion as a substantial
source of the observed discrepancy. Higher orders remain at physical contrast.
The calculation does not assign an additive percentage of the total error to Q,
or separately quantify the forward algorithm's omission of evanescent components.

## Numerical verification

| Check | Relative difference or residual |
|---|---:|
| Near vector cutoff, (109,55) vs (114,59) | 2.47e-12 |
| Near scalar cutoff | 2.72e-12 |
| Born potential derivative vs existing surface-integral function | 1.68e-05 |
| Born potential step, h=0.005 vs h=0.01 | 4.73e-05 |
| Second coefficient step refinement | 2.82e-05 |
| Batched Ex vs public spheroid_eval | 7.34e-16 |
| Scalar Cartesian gradient vs finite differences | 1.61e-08 |
| Exact logarithmic Helmholtz residual | 5.43e-07 |

Near-field cutoff checks use 257 independent grid points; public evaluation and
derivative/PDE checks use eight interior/exterior points at all four angles.
Born maps use a symmetric derivative in potential at h=0.005, verified against
`born_rytov_spheroid` with requested quadrature tolerance 1e-7. The detector
control uses a 51.2-um window with 0.1-um sampling, consistent with the earlier
sampling/window convergence audit in [oblate_report.md](oblate_report.md).
The figure was inspected after rendering; all panels use the same unclipped scale.

The run exposed a floating-point coordinate issue: four oblate axis points had
eta slightly larger than 1. `sph_coords.m` now clips eta_abs to its mathematical
upper bound before using it; the regression grid and existing coordinate and
vector-wave tests pass. No physical parameter or polarization convention changed.

## Files and reproduction

- `oblate_Ex_xz_4x3.png`, `.pdf`: requested figure.
- `oblate_xz_fields.mat`: complex Ex for all three columns, exact scalar field,
  incident field, Born field, Q, grid, geometry mask and cutoff diagnostics.
- `oblate_xz_exact.mat`: exact field and Cartesian scalar gradients.
- `oblate_xz_plot_errors.csv`: magnitude/complex errors over the whole window and interior.
- `oblate_rytov_cause.csv`: all 16 angle/contrast combinations.
- `oblate_xz_nonlinear_term.csv`, `oblate_xz_metrics.csv`: exact-field diagnostics.
- `oblate_xz_validation.csv`, `oblate_xz_derivative_validation.csv`: numerical checks.
- `oblate_xz.log`, `oblate_xz_derivative_validation.log`: successful execution evidence.

From the repository root in MATLAB:

```matlab
restoredefaultpath;
addpath('forward/exp');
oblate_na_sweep('xz');
oblate_na_sweep('xz_validate');
```

Then run `python forward/exp/plot_oblate_xz.py` in an environment containing
NumPy, SciPy and Matplotlib. This run used `/tmp/mie-oblate-plot-env/bin/python`.
MATLAB ran with `-nojvm -nodesktop -nosplash`; no global MATLAB settings changed.
