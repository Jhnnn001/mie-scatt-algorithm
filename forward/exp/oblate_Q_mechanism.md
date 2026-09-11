# Interpretation of the internal Q maximum

Read-only analysis of the existing exact fields on September 12, 2026.

The leading interpretation is concentration of internally reflected waves near
an axial caustic, followed by interference with the predominantly forward field.
The evidence supports this mechanism; it is not an exact decomposition by
internal-reflection order or a partition of the detector error.

At normal incidence the sampled maximum is at (x,z)=(0,-0.20) um:

- |U|=0.7041293921, normalized to unit scalar incident amplitude.
- grad(log|U/U0|)=(0,0,5.61290963) per um.
- grad(arg(U/U0))=(0,0,11.32853278) radians per um.
- The real amplitude-gradient, real phase-gradient, and imaginary cross terms
  divided by f are +2.82391705, -11.50331912, and +11.39902084 i.
- Q/f=-8.67940207+11.39902084 i; |Q|/f=14.32723618.

A Hann-windowed z Fourier transform over |z|<2.5 um at x=0 has its forward peak
at kz=16.1212 per um and its negative-kz peak at -9.22785 per um.
The negative/positive peak amplitude ratio is 0.11006; the fraction of squared
one-dimensional spectral norm at negative kz is 0.064984 on axis, 0.006553 at
x=0.32 um, and 0.000282 at x=2 um. These are finite-window spatial-spectrum
diagnostics, not optical reflected-power fractions. Applying the same window
to a pure forward plane wave gives a negative-kz squared-norm fraction of
2.22e-8. This supports a real backward-wave contribution concentrated near the axis.

Geometric rays obey Snell refraction at the lower entry boundary and specular
reflection at the upper boundary, with the actual nm=1.335381534 and np=1.365.
For normal incidence, representative reflected rays cross the z axis as follows:

| Entry radius (um) | Reflection x (um) | Reflection z (um) | Axis crossing z (um) |
|---:|---:|---:|---:|
| 1 | 0.984375 | 2.941286 | -1.038668 |
| 2 | 1.968731 | 2.757658 | -0.771842 |
| 3 | 2.953028 | 2.420881 | -0.240463 |
| 4 | 3.937078 | 1.849268 | 0.862068 |

See oblate_Q_reflected_ray_interpretation.png for these paths over the exact Q map.
The ray calculation includes neither Fresnel amplitude weights nor diffraction;
it explains the spatial concentration, not the exact field amplitude or Q peak.
The exact maximum can lie away from a geometric focus because Q is a logarithmic
field-gradient quantity and depends on interference with the forward field.

General background on internal-reflection caustics in dielectric scattering:
https://csuohio.elsevierpure.com/ws/portalfiles/portal/39955202/High-Order%20Interior%20Caustics%20Produced%20in%20Scattering%20of%20a%20Diagonal.pdf
That study is not a validation of the current geometry or parameter values.
