# Rytov detector-NA correction and spheroid comparison

Verified 2026-09-12 using MATLAB R2024a without the JVM.
The runnable checks are `test_born_rytov_voxel.m` and
`analyze_rytov_spheroid.m`; the numerical results and execution record are
`rytov_analysis.csv` and `rytov_analysis.log` in this directory.

## Corrected operator

The previous implementation restricted the linear Born field to detector NA
before using it to construct the Rytov phase. It also returned an exponentiated
field with nonzero frequencies outside that NA.

The corrected sequence is:

1. Sample the object's FFT on the forward propagating Ewald hemisphere.
2. Inverse-transform this unpupilled linear field and divide by the incident
   field to obtain the first Rytov complex phase.
3. Form the scattered Rytov field with `u0 .* expm1(phase)`.
4. Fourier-transform that field, apply detector NA, and inverse-transform.

Born retains its final linear pupil. Both returned scattered fields are now
bandlimited by the same detector pupil. `info.rytov_phase` describes the
pre-objective field: after pupil filtering it is generally different from
`log(1 + uR/u0)`.

The calculation includes forward propagating modes with `kz > 0`, omits
evanescent and exactly grazing modes, and retains the existing periodic
detector window and linear Ewald interpolation. Resolving the unpupilled
field requires `dx < lambda/(2*n_m)` and enough object-spectrum bandwidth
for the incident-shifted Ewald samples; these are now checked explicitly.
At the baseline wavelength and medium index, the first bound is 0.199194 um.
The baseline `dx=0.1 um` satisfies it; oblique illumination is additionally
checked against the available 3-D FFT coordinates.

## Physical parameters and comparison definition

| Quantity | Baseline |
|---|---:|
| Vacuum wavelength | 0.532 um (532 nm) |
| Wavelength in the medium | 0.398388 um |
| Medium / object index | 1.335381534 / 1.37, real |
| Axial / transverse semiaxis | a=5 um / b=2.5 um |
| Full dimensions | 10 x 5 x 5 um |
| Illumination / polarization | +z / lab x |
| Detector coordinate | z=7 um, 2 um beyond the forward tip |
| Detector NA | 0.10 |
| Voxel grid / spacing | 128 cubed / 0.1 um |
| Detector window | 12.8 x 12.8 um |
| FFT padding factor | 2 |
| Ray-estimated centerline phase delay | 4.08861 rad |

532 nm is an experimentally used live-cell ODT wavelength, not an arbitrary
unphysical scale: [Kim et al., 2014](https://doi.org/10.1364/OE.22.010398).
A 10-um-scale object is reasonable as a cell benchmark; it does not represent
all cell types, and this orientation deliberately puts its longest dimension
along the illumination. A cell flattened on a substrate can have a shorter
optical path. The sibling polarization experiment used the same baseline
indices, dimensions, and wavelength.

The pupil lies inside the sibling experiment's narrow angular square:
the inscribed circular limit is `n_m*sin(atan(0.08)) = 0.1064903`.
That experiment measured polarization leakage on a sampled patch; it did not
certify scalar amplitude/phase accuracy at every angle.

All reported errors below are complex **scattered-field relative L2 errors**:
`norm(model - exact)/norm(exact)`. They are not intensity errors, total-field
errors, or reconstruction errors. The earlier 187% / 33% values were for a
single spectral line; these tables use the full sampled circular pupil.

The reference is the validated vector solution from `spheroid_solve` and
`spheroid_eval`, projected onto Ex. For each detector spatial frequency,
the outgoing direction is `(kx,ky,kz)/k`. The far amplitude is extracted as
`A = r*exp(-i*k*r)*Ex_sca`, and converted to the detector spectrum using
`U = 2*pi*i*A*exp(i*kz*z_det)/kz`. Values at r and 2r, with
`r=1e6*max(a,b)`, agreed to at most 2.12e-6 relative error across the 15
comparisons. Each exact solution passed its `tol=1e-6` validation.

The independent continuous-shape control uses the closed-form spheroid
Fourier transform, `F(q)=f*V*3*(sin(Q)-Q*cos(Q))/Q^3`, where
`Q^2=b^2*(qx^2+qy^2)+a^2*qz^2` and the Q=0 limit is 1.
It removes voxelized-boundary and 3-D interpolation error, while preserving
the same first Rytov approximation and final NA. It is **not** an exact
scalar scattering solver. No call is made to `born_rytov_spheroid.m`.

## Results

| Baseline 128-cubed calculation | Relative error |
|---|---:|
| Born | 197.56% |
| Original returned Rytov field | 74.40% |
| Original Rytov field, with only a final pupil added | 26.21% |
| Corrected Rytov field: unpupilled phase, then final pupil | 25.74% |
| Continuous-shape first Rytov, same detector grid | 24.80% |

Most of the 74.40% to 25.74% reduction in this particular case comes from
removing out-of-pupil frequencies. Correcting the premature cutoff as well
changes the in-pupil error from 26.21% to 25.74%; the two changes must not be
confused. The premature cutoff can matter more for other windows: at a
25.6-um FOV, final-pupil-only gives 41.98%, versus 25.33% with the correct order.

| Numerical control, same physical object | Corrected voxel Rytov | Continuous-shape Rytov |
|---|---:|---:|
| Baseline: dx=0.1, padding=2, FOV=12.8 um | 25.74% | 24.80% |
| Padding=4 | 24.95% | 24.80% |
| dx=0.05, 256 cubed, same FOV and padding=2 | 25.78% | 24.80% |
| FOV=25.6 um, dx=0.1 | 25.33% | 24.00% |
| FOV=51.2 um, dx=0.05 | not run | 24.30% |
| FOV=102.4 um, dx=0.05 | not run | 24.34% |

The full complex-field difference between voxel and continuous-shape Rytov
is 3.41% at padding=2 and 1.44% at padding=4; these are distinct from the
difference between the two error percentages. Voxel/interpolation and finite
window errors cannot explain the persistent approximately 24% discrepancy.

| Physical control | Continuous-shape Rytov error |
|---|---:|
| Baseline, z=7 um, FOV=51.2 um | 24.30% |
| Same object, z=5.2 um | 16.98% |
| Same object, z=14 um | 47.34% |
| Same object, relative index contrast=0.001 | 1.04% |
| Same object, wavelength=633 nm | 22.84% |
| Half dimensions, z=7 um, FOV=12.8 um | 66.87% |
| Half dimensions, z=3.5 um, FOV=51.2 um | 24.21% |
| Half dimensions, z=2.7 um, FOV=51.2 um | 18.61% |

The wavelength control holds both indices fixed to isolate wavelength;
it is not a dispersive water/cell model at 633 nm. The half-size controls
show that reducing size while retaining the original detector coordinate
can worsen the result. Object size, wavelength, and propagation distance
jointly determine the diffraction pattern at the evaluation plane.

## Why first Rytov still differs

For `U = U0*exp(psi)` and plane-wave U0, the scalar Helmholtz equation becomes

`laplacian(psi) + 2*i*k_incident dot grad(psi) + grad(psi) dot grad(psi) = -f`.

First Rytov drops the quadratic gradient term. It can retain large accumulated
phase through the exponential without accurately describing strong spatial
variation of that complex phase. Smooth geometry alone does not make this
term small: this spheroid is thicker in the middle than at its edge, and its
boundary produces refraction and diffraction. Rytov's phase-gradient
condition and its advantages over Born for cell imaging are discussed in
[Sung et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2832333/).

An illustrative ray estimate gives the transverse phase profile
`phase(rho)=k0*(n_p-n_m)*2*a*sqrt(1-rho^2/b^2)`.
At `rho=b/sqrt(2)`, its slope is 1.635 rad/um, and the magnitude of its
squared transverse gradient divided by f is approximately 0.205. This is
only a scale estimate, not an exact validity test or a prediction of a
20.5% field error; it explains why an aspect-ratio-two homogeneous object
need not satisfy a small-phase-gradient approximation at this contrast.

There is also a direct numerical test of propagation consistency. Forming
first Rytov near the exit at z=5.2 um and then propagating the complex field
linearly through homogeneous water preserves its approximately 17% error
within the pupil. Evaluating the first Rytov formula anew at z=7 or 14 um
instead gives 24.30% or 47.34%. Direct and exit-propagated model fields differ
by 6.53% at z=7 and 26.93% at z=14. Exponentiating the propagated linear
phase is not equivalent to propagating the exponentiated field; the neglected
gradient term is what prevents exact free propagation of the first Rytov
field. This is a limitation of the approximation, separate from detector NA.
The near-exit control uses the same propagating-only model and does not
establish evanescent-field accuracy just 0.2 um from the surface.

These controls support index contrast and diffraction/propagation as major
contributors. They do not partition the residual into pure scalar-Rytov
error and scalar-versus-vector Maxwell error: that would require an exact
scalar spheroidal reference. Small discarded polarization alone cannot
justify attributing every percentage point to Rytov.

## Verification and compatibility

`test_born_rytov_voxel` first failed on the old implementation: changing
detector NA altered the retained Rytov field by 1.95%, and out-of-pupil
leakage was 0.778%. After the fix these values are 2.11e-16 and 1.86e-16.
Normal/tilted illumination, zero contrast, an independently summed translated
point source, and inadequate pre-objective bandwidth checks pass. The weak
spheroid comparison has 1.21% corrected voxel Rytov error and passes the 3%
regression bound. All 15 analysis rows completed; MATLAB Code Analyzer found
zero issues in the three changed/new MATLAB files, and `git diff --check`
passed. The analysis is a characterization, not a claim that all physical
cases meet an accuracy target.

The older `forward/tests/run_all.m` was not run because it calls the
explicitly excluded legacy spheroidal Born/Rytov model. Run the new checks
directly:

```matlab
addpath('forward/tests');
test_born_rytov_voxel;
results = analyze_rytov_spheroid;
```

The existing inverse test and default demo still use dx=0.2 um and the old
pre-exponential pupil convention. The inverse test now stops with
`born_rytov_voxel:RytovNyquistViolation`; its prior machine-precision equality
between physical Born and Rytov-derived data also needs revision because
logarithms and the physical pupil do not commute. These separately developed
inverse files were inspected but not edited or claimed to pass. Their
historical verification report predates this correction. Preserving the
12.8-um scene window with dx=0.1 requires 128 rather than 64 voxels per side.
