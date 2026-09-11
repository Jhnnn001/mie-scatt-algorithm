# Born and first Rytov for the flattened oblate case

Computed with MATLAB R2024a Update 9 on 2026-09-12.

The continuous-shape calculation gives an illumination-pupil-weighted mean
complex scattered-field error of **91.32% for Born and 10.51% for first Rytov**
against the co-polarized Maxwell reference. Rytov improves substantially over
the earlier prolate case, but this case does not meet a 5% accuracy target;
its weighted mean also slightly exceeds 10%. The scalar boundary-value model
itself has only **0.70%** mean co-polarized error.

## Conditions and polarization

| Quantity | Value |
|---|---:|
| Vacuum wavelength | 0.532 um |
| Background / object index | 1.335381534 / 1.365, real |
| Semiaxes | a=3 um along z, b=5 um along x and y |
| Full dimensions | 10 x 10 x 6 um |
| Detector plane | z=5 um, 2 um beyond the forward tip |
| Detector | Ideal circular pupil, NA=1.0, axis fixed along +z |
| Illumination NA | 0:0.05:0.5, theta=asin(NA/n_m), maximum 21.98877 degrees |
| Azimuth | 0:45:315 degrees; normal incidence counted once, 81 directions total |
| Original polarization | pol=[1,0] in incident_plane_wave.m |
| Baseline voxel calculation | 256 x 256 x 80, dx=0.1 um, FFT padding factor 2 |
| Baseline detector window | 25.6 x 25.6 um |

The retained convention is
`e0=Rz(phi) Ry(theta) Rz(-phi) xhat=cos(phi)e_TM-sin(phi)e_TE`.
For each theta, the batch calculation solves the xz-plane TM and TE systems;
at laboratory detector azimuth alpha it evaluates those solutions at
alpha-phi and combines them with cos(phi) and -sin(phi). The code checks that
the reconstructed incident vector matches the public incident_plane_wave
function. Every scalar/voxel incident direction is evaluated directly on the
same laboratory Cartesian FFT grid, including oblique azimuths.

## Error definitions and averaging

For each incident direction, `u_co=E_s dot conj(e0)` and
`epsilon_co=norm(u_model-u_co)/norm(u_co)` after the same detector pupil.
The norm covers the sampled complex scattered field; it is not an intensity,
total-field, phase-only, or reconstruction error. Parseval allows this norm
to be evaluated on the detector spectrum with the common frequency-cell area.

The mean assumes uniform area coverage of the illumination pupil, not equal
weight for each polar-angle ring. Ring edges are midpoints between the sampled
NA radii, truncated to [0,0.5]; each ring's weight is its annulus area divided
by the full pupil area, distributed equally among its azimuth samples.
The reported mean is the weighted mean of the per-direction errors. The
separately saved pooled L2 is the square root of the weighted total squared
error divided by the weighted total reference norm squared. No fields from
different illuminations are coherently averaged.

## Main results

All values below are percentages; maxima are over the 81 sampled directions.

| Calculation against co-polarized Maxwell reference | Weighted mean | Pooled L2 | Minimum | Maximum |
|---|---:|---:|---:|---:|
| Existing voxel Born, baseline grid | 89.679 | 89.716 | 88.095 | 91.968 |
| Existing voxel Rytov, baseline grid | 9.751 | 9.775 | 8.500 | 10.844 |
| Continuous-shape Born, refined 2-D grid | 91.319 | 91.368 | 89.084 | 93.648 |
| Continuous-shape Rytov, refined 2-D grid | 10.509 | 10.528 | 9.869 | 11.192 |
| Scalar boundary-value reference | 0.705 | 0.705 | 0.658 | 0.733 |

The refined continuous calculation uses the exact Fourier transform of the
homogeneous spheroid, a 51.2-um window and 0.05-um sampling. It retains the
same Born/Rytov equations as the voxel model, so it removes voxel boundary and
3-D Fourier-interpolation errors without removing the approximation error.
The Rytov phase is formed before the detector pupil, exponentiated, and then
filtered by that pupil, as in the corrected existing voxel implementation.

Against the scalar boundary-value reference, the continuous-shape mean errors
are **91.281% Born / 10.525% Rytov**. The near equality to the Maxwell
co-polarized comparison shows that the approximately 10.5% Rytov discrepancy
is predominantly in the scalar forward approximation, rather than in the
polarization reduction.

If each scalar result is converted to the vector `u_model e0` and compared
with the full Maxwell scattered field, the weighted means become
**91.327% Born / 10.935% Rytov** for continuous shape and
**89.689% Born / 10.210% Rytov** for the baseline voxel calculation.
These are scalar-model errors against a vector reference, not vector Born or
vector Rytov algorithms. The scalar boundary-value reference has **3.117%**
mean and **3.166%** maximum full-vector error, consistent with the earlier
oblate scalar-approximation result at this detector NA.

## Dependence on illumination angle

The following values average the available azimuths at each listed NA;
the full table also includes NA=0.05, 0.15, 0.25, 0.35 and 0.45.

| Illumination NA | Continuous Born co error (%) | Continuous Rytov co error (%) |
|---:|---:|---:|
| 0.0 | 89.084 | 9.869 |
| 0.1 | 89.252 | 9.921 |
| 0.2 | 89.766 | 10.065 |
| 0.3 | 90.646 | 10.310 |
| 0.4 | 91.923 | 10.680 |
| 0.5 | 93.648 | 11.187 |

The worst continuous Rytov co error is 11.192% at illumination NA=0.5,
phi=0 degrees. At fixed illumination NA, the maximum azimuthal spread of
continuous Rytov co error is only **0.0091 percentage points** in this case.
Azimuth was therefore a small numerical effect in the final error metric,
although invariance could not be assumed from the original polarization
definition alone. The voxel azimuth variation is larger because its Cartesian
grid and linear interpolation introduce additional direction dependence.

The ray-estimated normal-incidence centerline phase delay is 2.09885 rad,
compared with 4.08861 rad in the earlier prolate baseline. Reduced contrast
and optical thickness help, but the present object still accumulates a
substantial phase and Born remains inaccurate. The earlier reported corrected
voxel Rytov error was 25.74% for the prolate, normal-incidence, detector-NA=0.1
case; it is not a controlled comparison that isolates detector NA, because
index, shape, observation plane, detector aperture and illumination sampling
have all changed.

## Numerical audit

| Reference check | Relative difference or residual |
|---|---:|
| Vector cutoff: (109,55) versus (114,59) | 6.53e-13 |
| Scalar cutoff: (109,55) versus (114,59) | 6.75e-13 |
| Independent vector boundary nodes, maximum over incidence/polarization cases | 4.66e-12 |
| Independent scalar boundary nodes, maximum over incidence cases | 1.95e-12 |
| Analytic vector far limit versus public evaluator with radius extrapolation | 4.13e-8 |

Boundary residuals are weighted surface L2 residuals, not maximum pointwise
residuals. Vector checks include tangential E/H and normal n^2 E/H continuity;
unretained incident Fourier modes are included in the residual. Public
`spheroid_solve` runtime validation is not claimed for these experimental batch
solves. The small-size weak-contrast checks also compare both reference types
with analytic Born amplitudes, and the script checks the original polarization
convention and the absence of out-of-pupil voxel Rytov output.

For the continuous model, increasing the window from 25.6 to 51.2 um changes
the Rytov co-error metric by at most **0.0049 percentage points** across the
81 directions. On seven representative directions, direct spectral comparisons
at identical physical frequencies give a maximum **0.0739%** field difference
for 25.6 to 51.2 um, **0.0337%** for 51.2 to 102.4 um, and **0.000101%** for
dx=0.1 to 0.05 um at fixed 51.2-um window. These checks support an approximately
10.5% continuous Rytov error; extra printed digits document the calculation
and do not imply an accuracy bound at the solver's modal residual precision.

Coarsening illumination NA spacing from 0.05 to 0.1 changes the weighted mean
Rytov co error by **0.0110 percentage points**; using 90-degree instead of
45-degree azimuth spacing changes it by **0.000048 percentage points**.
These are sampling checks, not a rigorous bound over all continuous angles.

For the actual voxel algorithm, the numerical discrepancy from continuous
Rytov reaches **4.505% of the continuous Rytov field norm** at padding factor 2.
On the representative seven-direction control set, padding factor 3 reduces
the maximum discrepancy to **2.027%**, while changing the voxel co-error
metric by up to **0.748 percentage points** relative to padding factor 2.
Reducing voxel pitch to 0.0667 um at fixed window and padding factor 2 changes
that metric by at most **0.0363 percentage points**; increasing the detector
window to 51.2 um changes it by at most **0.0847 percentage points**.
Interpolation/padding error is therefore material in the baseline voxel result.
Its smaller 9.75% mean should not be interpreted as better physical accuracy
of the Rytov approximation: numerical and approximation errors partially cancel.
The voxel implementation has not been demonstrated converged below 1% field
error in this study.

## Scope and reproduction

The existing forward model uses only the forward propagating hemisphere;
evanescent and exactly grazing modes are omitted before forming the Rytov
phase. This audit quantifies the implemented model, and does not separately
isolate that omission from first-order Rytov truncation. The detector is an
ideal pupil with no objective polarization transport, aberration or analyzer;
actual camera cropping and measured illumination weights would define a
different experiment.

Run from the project root with MATLAB on the path:

```matlab
restoredefaultpath;
addpath('forward/exp');
oblate_na_sweep('check');
oblate_na_sweep('reference');
oblate_na_sweep('sweep');
oblate_na_sweep('fieldcheck');
oblate_na_sweep('report');
```

Numerical runs used `-nojvm -nodesktop -nosplash` with a session-local
`restoredefaultpath`; completed computation stages exited with status 0.
No global MATLAB configuration was changed.

- `oblate_na_sweep.m`: experiment, batched reference, independent checks and summaries.
- `oblate_summary.csv`: weighted means, pooled errors and extrema.
- `oblate_angles.csv`: all 264 baseline and convergence-control rows.
- `oblate_baseline.csv`, `oblate_window_51.csv`, `oblate_sampling_0p05.csv`: 81 directions each.
- `oblate_padding_3.csv`, `oblate_voxel_0p0667.csv`, `oblate_voxel_window_51.csv`: seven directions each.
- `oblate_reference_validation.csv`, `oblate_angular_convergence.csv`, `oblate_field_convergence.csv`: numerical audits.
- `oblate_reference_0.mat`, `oblate_reference_1.mat`: compact far-field reference coefficients and validation results; dense ODE trajectories are omitted after public-evaluator validation.
- `oblate_checks.log`, `oblate_reference.log`, `oblate_sweep.log`, `oblate_field_convergence.log`, `oblate_final_checks.log`: execution evidence.
