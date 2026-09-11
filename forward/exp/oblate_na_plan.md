# Oblate Born/Rytov evaluation

Approved scope: retain the existing `pol=[1,0]` convention; lambda=0.532 um,
n_m=1.335381534, n_p=1.365, axial semiaxis a=3 um, transverse b=5 um,
z_det=5 um, detector NA=1, illumination NA=0:0.05:0.5.

- [x] Add one reproducible experiment, `oblate_na_sweep.m`, with checkpointed
  reference and sweep stages; reuse the existing special functions and voxel
  forward model without changing production code.
- [x] Check the experimental batched scalar/vector boundary solver against
  weak-contrast analytic Born amplitudes and independent boundary nodes.
- [x] Solve TM and TE incidence in the xz plane for every polar angle at
  (L,M)=(109,55) and (114,59); verify cutoff convergence and the analytic
  far limit against the public `spheroid_eval` evaluator.
- [x] Compare continuous-shape and voxel Born/Rytov against scalar and vector
  references on a common detector spectrum. Recover the original x-polarization
  for azimuth phi by rotating cos(phi)*TM-sin(phi)*TE. Test direct azimuth
  evaluation for grid anisotropy. Apply the detector pupil after exponentiation.
- [x] Report per-angle co-polarized and full-vector relative L2 errors, scalar
  model error, uniform-illumination-pupil weighted mean, pooled L2 and worst case.
  Refine the 2-D window/sampling and voxel pitch/padding on representative angles.
- [x] Save CSV, MAT, execution logs, and a concise report with convergence
  limits. Do not interpret scalar-boundary or vector-reference agreement as
  Rytov accuracy, and do not claim public runtime validation for batched solves.

Run with MATLAB R2024a, no JVM, a session-local `restoredefaultpath`, and
`addpath('forward/exp'); oblate_na_sweep('check')`, then `'reference'`,
`'sweep'`, `'fieldcheck'`, and `'report'`. The script's assertions are the runnable checks; no test framework
or dependency is added. This is an in-place scientific experiment; no commits,
worktree changes, or modifications to earlier results are required.
