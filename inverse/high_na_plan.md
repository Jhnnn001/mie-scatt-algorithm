# Matched high-NA ODT implementation and experiment plan

User approved on 2026-09-13: execute in this existing `inverse/` folder through
the saved numerical conclusion. Do not change the other section or Git history.

**Goal:** Compare zero-filled Fourier, nonnegative Gerchberg-Papoulis (GP),
and Lim TV reconstructions at illumination NA 0.5 and detector NA 1.0.

**Architecture:** Retain nonuniform Ewald samples. Use a padded unitary FFT
and trilinear sampling with its exact transpose. Physical detector fields come
from the existing validated oblate reference helpers, independently of the
voxel inverse. A matched linear control separates numerical and model error.

**Tools:** Existing MATLAB R2024a and NumPy/SciPy/Matplotlib; no installs.

## Conditions

- Lambda 0.532 um, water index 1.335381534; a=3 um (z), b=5 um (x/y).
- Latest user selection: only t=0.01, n_p=1.33567771866,
  delta_n=0.00029618466; +z detector at z=5 um. The other contrast cases
  remain forward-study provenance, not additional inverse experiments.
- Retain the 81 existing directions, pol=[1,0], and Maxwell co-polarization.
- Initial detector 512x512 and volume 128x128x80, pitch 0.1 um.
- Real nonnegative chi=n^2/n_m^2-1; no true shape/support in reconstruction.
- No added noise. Lower contrast also means weaker signal, not fixed-noise SNR.

## Execution checklist

- [x] **Matched operator/preparation:** create `odt_operator.m`, update
  `odt_prepare.m` for independent detector/volume sizes, post-pupil logarithms,
  and raw samples. Test complex adjoint identity <1e-11 and agreement with
  the existing voxel forward <1e-11; retain analytic point-source/phase tests.
- [x] **Three solvers:** update `odt_reconstruct.m`. Direct uses normalized
  adjoint gridding, Hermitian completion and zero fill. GP alternates a
  minimum-norm least-squares data correction (CGLS) and positivity. TV uses
  Lim Eq.14 inner splitting and Eq.15 residual reinjection, with the true AHA
  system solved by CG. One outer iteration retains penalized TV as a control.
  Tests: full-data identity, exact two-level TV, outer reinjection, zero data,
  capped-iteration status and matched incomplete-data behavior.
- [x] **Demonstration:** update `demo_inverse.m` to valid sampling and the
  oblate/index defaults. Physical Rytov data are not an exactly linear model.
  Run the small example through `inverse/tests/run_all.m`.
- [x] **Experiment:** create `run_inverse_experiment.m`; reuse
  `oblate_na_sweep('helpers')` and completed contrast reference checkpoints.
  Check exact parameter identity and reference validation. Save sample and
  reconstruction checkpoints, settings, CSV metrics and histories in
  `inverse/results/`. Use the same exact fields for all three methods. Include
  a matched linear control, padding/FOV/sampling and iteration/TV sensitivity.
  Primary regularization uses measured data scale, not the true object.
- [x] **Analysis/review:** save figures and Korean `results/report.md`, update
  `README.md` and `verification.md`. Report contrast-normalized RI error,
  mean RI bias, axial half-contrast span, background error, Ewald/physical-field
  residuals, convergence and time. Run tests, Code Analyzer, CSV/figure checks,
  independent review and affected regressions before completion.

## Reproduction

```matlab
restoredefaultpath
run('inverse/tests/run_all.m')
addpath('inverse')
settings = struct('padding',4,'gp_iter',100,'tv_inner',200,'cg_max_iter',100);
run_inverse_experiment('exact',settings)
```

No commit or push is part of this experiment.

## Final run inventory

Finish and analyze these 17 reconstructions; do not expand the parameter search
indefinitely if a cap or model limitation remains. Numerical nonconvergence is
reported as a result, never relabeled as a successful solve.

| Label | Methods | Diagnostic |
|---|---|---|
| exact | direct, gp, tv | Independent physical data, padding 4 |
| linear_control | direct, gp, tv | Matched A*truth, padding 4 |
| padding8 | direct, gp, tv | Fourier interpolation refinement |
| gp_control_tight | gp | Matched data, 300 cap, CGLS tolerance 1e-5 |
| gp_tight | gp | Exact data, same tighter GP settings |
| tv_long | tv | 1000 inner cap, baseline TV weight and penalties |
| tv_half | tv | Half the TV weight, longer inner solve |
| tv_double | tv | Double the TV weight, longer inner solve |
| tv_penalty | tv | Same objective, rho=rho_positive=0.01 |
| tv_refined | tv | Padding 8, physically scaled splitting penalties |
| tv_strict | tv | Same refined TV problem, tolerance 1e-5 |

The penalty checks were added after the longer run exposed significant
finite-tolerance sensitivity. They change numerical convergence, not the TV
objective. The strict check is the final numerical sensitivity experiment.

Completed 2026-09-13: all 11 runs / 17 reconstructions, sampling checks and
11 figures are saved. Seven iterative runs remain explicitly unconverged;
completion means the declared experiments finished, not that every solver
converged or accurate RI recovery was achieved. See `results/artifact_check.log`
and `results/verification_final.log` for the final checks.
