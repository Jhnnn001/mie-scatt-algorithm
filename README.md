# Forward scattering and ODT reconstruction

MATLAB experiments comparing Born/Rytov approximations with a full-vector Maxwell spheroid solver and reconstructing 3D refractive index by direct Fourier inversion, Gerchberg–Papoulis (GP), and total variation (TV). GP and TV impose nonnegative index contrast.

## Dependencies

- MATLAB R2024a; base MATLAB functions, no additional toolboxes or Java required for the commands below.
- Python 3.11, NumPy 2.4, SciPy 1.17, and Matplotlib 3.11 for checking saved results and drawing figures.
- Conda for creating the Python environment.

## Setup

```sh
git clone https://github.com/Jhnnn001/mie-scatt-algorithm.git
cd mie-scatt-algorithm
conda env create -f environment.yml
conda activate mie-odt
```

Make `matlab` available on your shell's PATH; for the macOS R2024a installation:

```sh
export PATH="/Applications/MATLAB_R2024a.app/bin:$PATH"
```

If MATLAB cannot find built-in functions such as `addpath`, restore its factory search path for this shell session:

```sh
export MATLABPATH="$(perl /Applications/MATLAB_R2024a.app/toolbox/local/getphlpaths.pl /Applications/MATLAB_R2024a.app)"
```

## Run

Check the solvers and regenerate the four figures from the included results:

```sh
matlab -nojvm -batch "run_project('check')"
python plot_results.py
```

Optional literature benchmark:

```sh
matlab -nojvm -batch "addpath('maxwell-solver','maxwell-solver/tests'); test_barton2001"
```

Reproduce the experiments:

```sh
matlab -nojvm -batch "run_project('forward')"  # Prolate baseline, oblate size/contrast sweeps
matlab -nojvm -batch "run_project('fields')"   # Exact/Born/Rytov sections and the omitted Q term
matlab -nojvm -batch "run_project('inverse')"  # Direct, GP, TV, and matched-data control
python plot_results.py
```

Saved reconstructions are reused. Delete `reconstruction/results/<run>/<method>.mat` to recompute one; figures are saved in `figures/`.

## Project structure

```text
maxwell-solver/             Full-vector spheroid solver and validation tests
born-rytov/                 Born/Rytov forward models and tests
experiments/
  prolate-baseline/         Prolate forward comparison and error table
  size-sweep/               Semiaxis scales of 100%, 10%, and 1%
  contrast-sweep/           Index-contrast scales of 100%, 50%, and 1%
  field-analysis/           Exact/Born/Rytov fields and the omitted Q term
  run_parameter_sweep.m     Shared size/contrast experiment runner
reconstruction/            Ewald sampling, inverse solvers, and tests
  results/
    maxwell-direct-gp/      Direct and GP results from Maxwell data
    maxwell-tv/             TV result from Maxwell data
    matched-control/       Direct/GP/TV results from the discrete forward model
figures/                   Four generated result figures
```

## Results

| Experiment | Relative error | Reference |
|---|---|---|
| Prolate baseline | Born **197.6%**, Rytov **25.7%** | Maxwell `Ex` |
| Oblate baseline | Born **91.3%**, Rytov **10.5%** | Co-polarized Maxwell field |
| Size reduced to 1/100 | Born and Rytov **0.72%** | Scalar Helmholtz field |
| Index contrast reduced to 1/100 | Born **0.88%**, Rytov **0.091%** | Scalar Helmholtz field |
| Low-contrast oblate reconstruction | Direct **46.9%**, GP **36.4%**, TV **28.9%** | True RI volume |

Forward errors are relative L₂ errors of the complex scattered field within the detector pupil; oblate values are illumination-pupil-weighted means over 81 directions.
Reconstruction error is `||n_reconstructed − n_true||₂ / ||n_true − n_background||₂` over the full volume.
At normal incidence, the scalar Rytov term `Q = ∇ψ · ∇ψ` reaches `|Q|/|f| = 14.33` inside the oblate xz slice, showing that the omitted term is locally large ([analysis](experiments/field-analysis/rytov_nonlinear_term.md)).

![Direct, GP, and TV reconstruction](figures/reconstruction.png)

White contours mark the true boundary.

The reported cases use vacuum wavelength 0.532 μm and background index 1.335381534, with polarization defined by `pol=[1,0]` in `incident_plane_wave.m`.

| Case | Semiaxes (z, x/y), μm | Object index | Max illumination NA / detector NA | Detector z, μm |
|---|---|---|---|---|
| Prolate | 5, 2.5 | 1.37 | 0 / 0.1 | 7 |
| Oblate | 3, 5 | 1.365 | 0.5 / 1.0 | 5 |
| Inverse | 3, 5 | 1.33567771866 | 0.5 / 1.0 | 5 |

The inverse uses a 128 × 128 × 80 volume at 0.1 μm spacing; its low-contrast Rytov forward error against the co-polarized Maxwell field is 0.66%.
The reported GP and TV results have weighted data residuals below 1% despite substantial RI error; missing-cone ambiguity, discretization, and regularization can all contribute.
GP is unconverged after 100 iterations. These are homogeneous-spheroid simulations with plane-wave illumination, an ideal pupil, and no added noise.

## References

- J. P. Barton, [*Applied Optics* **40**, 3598–3607 (2001)](https://doi.org/10.1364/AO.40.003598): spheroidal Maxwell fields.
- S. Asano and G. Yamamoto, [*Applied Optics* **14**, 29–49 (1975)](https://doi.org/10.1364/AO.14.000029): spheroidal wave-function expansion.
- A. C. Kak and M. Slaney, [*Principles of Computerized Tomographic Imaging*, Chapter 6 (1988)](https://www.slaney.org/pct/): diffraction tomography.
- J. Lim et al., [*Optics Express* **23**, 16933–16948 (2015)](https://doi.org/10.1364/OE.23.016933): missing-cone regularization.
- T. Goldstein and S. Osher, [*SIAM Journal on Imaging Sciences* **2**, 323–343 (2009)](https://doi.org/10.1137/080725891): split Bregman optimization.
