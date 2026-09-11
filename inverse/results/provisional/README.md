# Provisional diagnostics, not final reconstruction results

These early runs used implementation 2: the direct Fourier gridding grid
changed with internal interpolation padding. This also changed initialization
and the data-derived TV parameter. The corrected implementation uses a fixed
direct reconstruction grid and preserves the physical TV/data balance across
padding and voxel-pitch refinements.

The padding2 `baseline` completed and identified substantial interpolation
error. The old `matched` and `padding4` runs were intentionally terminated
after that gridding issue was reproduced. Their partial checkpoints/logs are
preserved for provenance and must not be mistaken for completed final runs.

Direct/GP stationarity entries of zero in these old files were placeholders,
not measured KKT residuals. The current implementation reports these as NaN.
