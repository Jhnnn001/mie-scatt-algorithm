# Index-contrast experiment plan

Goal: measure whether reducing the index contrast makes first Rytov accurate, with a separate reproducible dataset and Korean readme.md.

- [x] Reuse ../1/parameter_sweep.m and the existing oblate solver. Keep a=3 um, b=5 um and all optical/detector parameters from experiment 1.
- [x] Set n_p = n_m + t*(1.365-n_m), with t = 1, 0.5, 0.01. This scales delta n, not n_p squared minus n_m squared. The smallest ratio remains strictly greater than one.
- [x] Evaluate scalar Rytov error against scalar exact and Maxwell references on 81 illumination directions. Retain Born and vector-Born controls so truncation and scalar-model errors can be separated.
- [x] Verify references and spatial/angular convergence; save CSV, MAT, execution logs and plots in this folder. Reuse the common current-condition checkpoint where its parameters match.
- [x] Write readme.md with exact index values/ratios, phase scales, the 1% and 5% error criteria, computed results, limitations and reproduction instructions. Finish all checks before reporting completion.

Implementation: run_contrast_experiment.m calls ../1/parameter_sweep.m with project='contrast'. Plotting/report data use the existing Python environment. No production forward or inverse model is replaced.
