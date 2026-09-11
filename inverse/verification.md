# Verification record

## Current high-NA implementation, 2026-09-13

The historical section below describes the replaced nearest-bin implementation,
not the current experiment. Current settings and results are in `README.md`
and `results/report.md`.

- The full inverse suite passed on the final implementation, including the
  exact separable Fourier optimization and runner penalty settings
  (`results/verification_final.log`):
  physical pupil/log handling, analytic FFT signs/scales, complex adjoint,
  full-FFT equivalence, direct-grid padding invariance, analytic TV endpoints,
  Lim reinjection, independent real-pseudoinverse GP, stationarity false-stop
  regression, exact oblique detector energy, and existing voxel-forward match.
- MATLAB Code Analyzer reported 0 issues across all 12 inverse MATLAB files
  in that same final run.
- The accepted-Nyquist-roundoff crop regression was reproduced as an
  out-of-range sparse index, then fixed by using the interpolation's same
  endpoint clamping. The full inverse suite passed afterward; the separate
  objective/penalty scaling-invariance regression also passed.
- Numerical sampling checks completed: padding 2/3/4/6/8, voxel pitch
  0.1/0.08 μm, and independent detector pitch/window changes on 7 angles.
- A separate run of `forward/tests/run_all.m` stopped at
  `test_born_rytov.m:178`: the test expects
  `born_rytov_spheroid:InvalidTolerance` for `tol=99*eps`, but the current
  function validates only positive tolerance and has no such rejection.
  This unrelated forward function/test mismatch was not edited; therefore
  a repository-wide passing test suite is **not** claimed.
- A read-only independent review found and prompted fixes for hidden TV
  stationarity error, insufficient GP testing, cache/label validation, and
  exact-data provenance checks; its follow-up found no substantive defect in
  the reduced Fourier operator, native direct normalization, or physical TV
  scaling. CGLS projection accuracy and capped ADMM iterations remain explicit
  numerical limitations requiring the saved sensitivity comparisons.
- A final read-only audit matched all 11 run configurations to the reproduction
  commands and checked the completed numerical tables and new plot branches.
  The strict-TV and sensitivity figures were then generated and visually
  inspected; final independent metric recomputation passed for all 17
  reconstructions. `results/artifact_check.log` confirms all 11 run inventories,
  matching MAT/CSV settings, histories and convergence flags, 11 figures and
  21 local document links after final report edits. Seven iterative runs remain
  explicitly unconverged.
- `results/source_data.sha256` records the 22 inverse source/reference/cache
  files used for provenance; all recorded hashes verified successfully.

Some diaries contain interrupted attempts from mathematically equivalent FFT
optimizations. Completed MAT checkpoints and their CSV histories identify the
reported reconstructions; the final PASS marker identifies run completion,
not a successful numerical convergence flag for every method.

## Historical verification, 2026-09-12

MATLAB R2024a Update 9 on macOS, `-nojvm`, using the process-local search-path
setup in README.md. No existing forward/exact solver files were edited.

## Numerical checks

Final rerun: `tests/run_all.m` passed. MATLAB Code Analyzer reported zero
issues across all five MATLAB files in `inverse/` and `inverse/tests/`.

`tests/run_all.m` checks:

- Central voxel: recovered continuous Fourier samples agree with chi0*dx^3
  within 2e-12, including tilted illumination.
- Analytic displaced point scatterer: recovered q coordinates and complex
  samples agree with the independently synthesized spectrum within 1e-12
  and 1e-11 respectively.
- Born/Rytov sample agreement, arbitrary nonzero complex reference amplitude,
  changing detector NA, the DC-only low-NA limit and invalid reference rejection.
- A smooth 4.5-radian phase bump is recovered automatically and with supplied
  phase; inconsistent supplied phase and unsupported border phase are rejected.
- Complete-data inverse transform agrees to relative 1e-12; GP preserves an
  already nonnegative complete-data volume.
- Periodic two-level TV solution agrees with the analytic endpoints 0.034 and
  0.086 within 2e-6 (input 0.03/0.09, alpha=0.008, four voxels per level).
- Zero signal and alpha=0 remain finite; a one-iteration TV run correctly
  reports unconverged; a water/cell integration case exercises all three paths.

## Demonstration results

All use wavelength 0.532 um, water 1.335381534, cell 1.37, no noise.
Relative error is ||n_reconstructed-n_true|| / ||n_true-n_water|| over the
whole volume; residual is on original Fourier samples, before nearest-bin
gridding. `converged` is a numerical criterion, not scalar validity.

Small integration case: 32^3, dx=0.2 um, a=1 um, b=0.6 um, 13 angles,
NA_illum=0.2, NA_det=0.55, Rytov, alpha=1e-4, max_iter=100.

| Method | Relative RI error | Raw sample residual | Converged |
|---|---:|---:|---|
| Direct | 0.74819 | 0.096145 | Yes (direct transform) |
| GP | 0.40811 | 0.036563 | No |
| TV | 0.30358 | 0.015663 | No |

Full-size case: 64^3, dx=0.2 um, a=5 um, b=2.5 um, 49 angles,
alpha=1e-4. These cases also show why improvement must be measured rather
than assumed.

| Model | NA_det | NA_illum | Method | Relative RI error | Raw sample residual | Converged / cap |
|---|---:|---:|---|---:|---:|---|
| Rytov | 0.7 | 0.2 | Direct | 0.60708 | 0.099727 | Direct |
| Rytov | 0.7 | 0.2 | GP | 0.60196 | 0.056473 | No / 500 |
| Rytov | 0.7 | 0.2 | TV | 0.63434 | 0.054537 | No / 500 |
| Born | 0.4 | 0.2 | Direct | 0.62402 | 0.098853 | Direct |
| Born | 0.4 | 0.2 | GP | 0.55287 | 0.056012 | No / 500 |
| Born | 0.4 | 0.2 | TV | 0.68692 | 0.053365 | No / 500 |
| Rytov | 0.7 | 0.5 | Direct | 0.46433 | 0.10280 | Direct |
| Rytov | 0.7 | 0.5 | GP | 0.78062 | 0.072547 | No / 3000 |
| Rytov | 0.7 | 0.5 | TV | 0.51992 | 0.070875 | Yes, 745 iterations |

The last TV case stopped at primal residual 2.97e-6 and dual residual
9.98e-6. Its RI error still exceeds the direct reconstruction's error;
numerical convergence is not an accuracy guarantee. Nearest-bin gridding,
finite FOV and limited angular support are not removed by running longer.
The 3000-step GP case ended with relative iterate change 9.56e-5, above
tol=1e-5. No scalar-valid NA or successful experimental reconstruction is
inferred from these simulations.

An independent read-only review checked Fourier normalization, Hermitian
indexing, TV gradient/adjoint signs, shrinkage and positivity updates, the
stopping diagnostics, and demo metrics. The displaced-scatterer and wrapped-
phase tests were added in response to that review. Plotting was not exercised in the verification session.
