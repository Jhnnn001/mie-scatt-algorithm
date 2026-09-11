# Size experiment plan

Goal: measure Born and Rytov scattered-field error for three oblate sizes, with reproducible data and a Korean readme.md.

Authorized by the user on 2026-09-13: finish both experiments in folders 1 and 2. Execute in this workspace; do not commit or change earlier results.

- [x] Reuse the validated batched reference and continuous Born/Rytov formulas through a minimal helper interface in forward/exp/oblate_na_sweep.m. Check the weak-contrast limit and polarization before the sweep.
- [x] Keep wavelength 0.532 um, background index 1.335381534, particle index 1.365, a:b = 3:5, z_det = 5 um, illumination NA = 0:0.05:0.5, detector NA = 1. Scale both semiaxes by 1, 0.1, 0.01.
- [x] Solve converged scalar and Maxwell references for every size. Compare scalar Born/Rytov against scalar exact, Maxwell co-polarized and full-vector fields; include transverse-projected vector Born as a control for scalar-model breakdown.
- [x] Use 81 illumination directions and uniform illumination-pupil area weights. Save per-angle CSV, weighted/pooled/worst summaries, reference checkpoints, and logs.
- [x] Check modal cutoff, independent boundary nodes, the public evaluator, detector window/sampling, and angular coarsening. Keep convergence limits visible in the report.
- [x] Plot size dependence and angle dependence; write readme.md with parameters, formulas, results, interpretation, validation and reproduction commands. Verify CSV dimensions, plots and links.

Implementation: run_size_experiment.m calls parameter_sweep.m in this folder. The contrast experiment reuses this runner. Expose existing local functions rather than copying their solver implementation. Checks are executable MATLAB assertions plus independent CSV/plot verification; no test framework or new dependency.
