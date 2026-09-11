# Linear ODT Inverse Implementation Plan

Historical initial implementation; superseded by `high_na_plan.md` on 2026-09-13.

> **For agentic workers:** Implement sequentially in this session; use the
> test-first and verification workflows. All task files belong in `inverse/`.

**Goal:** Reconstruct a real, isotropic, nonabsorbing 3-D refractive-index
volume from Born or Rytov complex detector fields with variable detector NA.

**Architecture:** Reuse `forward/born_rytov_voxel.m` for simulated measurements
and geometry validation. Convert total fields to Ewald samples, average onto
a Cartesian Fourier grid, and run zero fill, Gerchberg-Papoulis (GP), or
nonnegative isotropic TV with Split Bregman on that same gridded data.

**Tech Stack:** Base MATLAB R2024a, no additional dependencies.

## Constraints and decisions

- Work in the requested `inverse/` folder without changing existing solvers.
- Vacuum wavelength 0.532 um; water index 1.335381534; cell index 1.37.
- Illumination travels with positive kz; fixed detector is beyond the +z grid face.
- Detector NA is an input, not a certified scalar-validity threshold.
- Illumination angles stay independent of detector NA, and must fit its pupil.
- Unknown is chi = n^2/n_m^2 - 1; TV and positivity act on chi, not n itself.
- Use unitary centered 3-D FFTs and periodic, unit-voxel forward differences.
- Nearest-bin Ewald gridding is a documented approximation; validate Fourier
  scaling separately and expose raw sample locations and gridding displacement.
- Real chi implies Hermitian Fourier symmetry; conjugate bins are inferred
  from that assumption, not additional measured angles.
- Rytov input is total field, not scattered field. Automatic 2-D unwrapping
  checks both axis orders and assumes a near-zero background border phase;
  caller-supplied unwrapped phase is supported and checked modulo 2*pi.
- The current forward's pupil acts before the Rytov exponential. Reusing it
  verifies that effective linear model, not physical post-field pupil filtering.
- TV solves 0.5*||M F chi-y||_2^2 + alpha*TV(chi), chi >= 0.
  This uses Goldstein-Osher splitting and Lim Appendix A's gradient/positivity
  updates without Lim Eq. (15)'s outer exact-data Bregman reinjection.
- Use no true support or true shape in reconstruction. Report incomplete
  numerical convergence instead of treating an iteration limit as success.

## Tasks and verification

1. **Data preparation** — `odt_prepare.m` consumes total U, reference u0,
   forward geometry info, and 'born'/'rytov'; returns gridded unitary spectrum,
   mask, sampling diagnostics and geometry. Test an isolated central voxel:
   its exact Fourier samples equal chi0*dx^3 at every q; this catches propagation
   sign, reference normalization, and FFT scaling errors. Test varying NA,
   off-axis illumination, Born/Rytov agreement and zero reference rejection.
2. **Three solvers** — `odt_reconstruct.m` consumes prepared data and method
   'direct'/'gp'/'tv'; returns n, chi and solver diagnostics. Test complete-data
   identity; GP nonnegativity; and an independently solvable periodic two-level
   TV problem (four voxels per level: endpoints shift by alpha/2).
3. **Runnable example** — `demo_inverse.m` reuses the existing voxel forward,
   supports separate NA_det/NA_illum and both models, reports error in n-n_m,
   mean object RI, axial extent, raw linear-field residual and solver status.
   Run all three methods for water/cell data and more than one NA.
4. **Documentation and review** — `README.md` documents equations, conventions,
   polarization scope, command examples, references and measured validation.
   Run `tests/run_all.m`, MATLAB `checkcode`, and an independent numerical/code
   review; fix actionable findings and rerun affected checks.

Test command:

```sh
/Applications/MATLAB_R2024a.app/bin/matlab -nojvm -batch "run('inverse/tests/run_all.m')"
```

## Status

- [x] Failing tests observed before implementation.
- [x] Data preparation checks pass.
- [x] Solver checks pass.
- [x] Water/cell demonstration and variable-NA checks pass (finite outputs;
      reconstruction accuracy and numerical convergence are reported separately).
- [x] Independent numerical/code review completed; final verification recorded
      in `verification.md`.
