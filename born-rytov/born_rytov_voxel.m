function [uB, uR, u0, info] = born_rytov_voxel( ...
    lambda, n, n_m, dx, theta_inc, phi_inc, NA_det, z_det, opts)
%BORN_RYTOV_VOXEL Scalar FFT/Ewald Born and Rytov detector fields.
% n(ix,iy,iz) is sampled on an even, centered Cartesian grid with isotropic
% voxel pitch dx. Light travels toward +z and the detector plane is z=z_det.
% theta_inc and phi_inc are radians and may be vectors or scalar-expanded.
% uB and uR are scattered fields of size Nx-by-Ny-by-Nillum; u0 is incident.
% The detector has the same lateral pitch and field of view as the voxel grid.
% Both outputs include the detector pupil. Rytov is exponentiated BEFORE
% that pupil: info.rytov_phase is the pre-objective phase, not log(1+uR/u0).
% The pre-objective field uses the forward propagating hemisphere (kz>0);
% evanescent and exactly grazing modes are omitted. dx must resolve this
% hemisphere and its incident-shifted Ewald samples, even for a small NA.

started = tic;
if nargin < 9
    opts = struct;
end

validateattributes(lambda, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename, 'lambda');
validateattributes(n, {'numeric'}, {'nonempty','finite'}, mfilename, 'n');
validateattributes(n_m, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename, 'n_m');
validateattributes(dx, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename, 'dx');
validateattributes(theta_inc, {'numeric'}, ...
    {'real','finite','nonempty'}, mfilename, 'theta_inc');
validateattributes(phi_inc, {'numeric'}, ...
    {'real','finite','nonempty'}, mfilename, 'phi_inc');
validateattributes(NA_det, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename, 'NA_det');
validateattributes(z_det, {'numeric'}, ...
    {'real','finite','scalar'}, mfilename, 'z_det');
if ndims(n) ~= 3
    error('born_rytov_voxel:InvalidVolume', 'n must be a three-dimensional array.');
end
grid_size = [size(n,1), size(n,2), size(n,3)];
if any(grid_size < 4) || any(mod(grid_size,2) ~= 0)
    error('born_rytov_voxel:InvalidGridSize', ...
        'Every dimension of n must be even and at least 4.');
end
if ~isstruct(opts) || ~isscalar(opts)
    error('born_rytov_voxel:InvalidOptions', 'opts must be a scalar struct.');
end
known = {'pad_factor'};
if any(~ismember(fieldnames(opts), known))
    error('born_rytov_voxel:UnknownOption', 'opts contains an unknown field.');
end
if ~isfield(opts, 'pad_factor')
    opts.pad_factor = 2;
end
validateattributes(opts.pad_factor, {'numeric'}, ...
    {'real','finite','scalar','integer','positive'}, mfilename, 'opts.pad_factor');

lambda = double(lambda);
n = double(n);
n_m = double(n_m);
dx = double(dx);
NA_det = double(NA_det);
z_det = double(z_det);
if NA_det >= n_m
    error('born_rytov_voxel:InvalidNA', ...
        'NA_det must be smaller than n_m so every collected mode has positive kz.');
end

theta_inc = double(theta_inc(:));
phi_inc = double(phi_inc(:));
if isscalar(theta_inc)
    theta_inc = repmat(theta_inc, size(phi_inc));
elseif isscalar(phi_inc)
    phi_inc = repmat(phi_inc, size(theta_inc));
elseif numel(theta_inc) ~= numel(phi_inc)
    error('born_rytov_voxel:IncidentAngleSizeMismatch', ...
        'theta_inc and phi_inc must have equal lengths or one must be scalar.');
end
if any(theta_inc < 0 | theta_inc >= pi/2)
    error('born_rytov_voxel:InvalidIncidentDirection', ...
        'theta_inc must satisfy 0 <= theta_inc < pi/2 for +z illumination.');
end

k0 = 2*pi/lambda;
k = n_m*k0;
cutoff = k0*NA_det;
if cutoff >= pi/dx
    error('born_rytov_voxel:DetectorNyquistViolation', ...
        'dx must satisfy k0*NA_det < pi/dx for the requested detector NA.');
end
if k >= pi/dx
    error('born_rytov_voxel:RytovNyquistViolation', ...
        'Rytov before the pupil requires dx < lambda/(2*n_m).');
end
k_incident = k * [sin(theta_inc).*cos(phi_inc), ...
    sin(theta_inc).*sin(phi_inc), cos(theta_inc)];
transverse_incident = hypot(k_incident(:,1), k_incident(:,2));
if any(transverse_incident > cutoff + 64*eps(k))
    error('born_rytov_voxel:ReferenceOutsidePupil', ...
        'Every incident reference wave must lie inside the detector NA.');
end

x = (-grid_size(1)/2:grid_size(1)/2-1) * dx;
y = (-grid_size(2)/2:grid_size(2)/2-1) * dx;
z = (-grid_size(3)/2:grid_size(3)/2-1) * dx;
if z_det <= z(end)+dx/2
    error('born_rytov_voxel:DetectorInsideGrid', ...
        'z_det must be above the +z face of the voxel grid.');
end

f = k0^2 * (n.^2 - n_m^2);
padded_size = double(opts.pad_factor) * grid_size;
f_padded = zeros(padded_size, 'like', f);
offset = (padded_size-grid_size)/2;
ix = offset(1) + (1:grid_size(1));
iy = offset(2) + (1:grid_size(2));
iz = offset(3) + (1:grid_size(3));
f_padded(ix,iy,iz) = f;
F = fftshift(fftn(ifftshift(f_padded))) * dx^3;
clear f_padded

dk_padded = 2*pi ./ (padded_size*dx);
kx_padded = (-padded_size(1)/2:padded_size(1)/2-1) * dk_padded(1);
ky_padded = (-padded_size(2)/2:padded_size(2)/2-1) * dk_padded(2);
kz_padded = (-padded_size(3)/2:padded_size(3)/2-1) * dk_padded(3);
kx = (-grid_size(1)/2:grid_size(1)/2-1) * (2*pi/(grid_size(1)*dx));
ky = (-grid_size(2)/2:grid_size(2)/2-1) * (2*pi/(grid_size(2)*dx));
[KX, KY] = ndgrid(kx, ky);
transverse_squared = KX.^2 + KY.^2;
pupil = transverse_squared <= cutoff^2;
propagating = transverse_squared < k^2;
KZ = nan(size(KX));
KZ(propagating) = sqrt(k^2-transverse_squared(propagating));
propagator = complex(zeros(size(KX)));
propagator(propagating) = 1i * exp(1i*KZ(propagating)*z_det) ./ ...
    (2*KZ(propagating));

n_illumination = numel(theta_inc);
output_size = [grid_size(1), grid_size(2), n_illumination];
uB = complex(zeros(output_size));
uR = complex(zeros(output_size));
u0 = complex(zeros(output_size));
phiR = complex(zeros(output_size));
born_spectrum = complex(zeros(output_size));
[X, Y] = ndgrid(x, y);
interpolation_tolerance = 64*eps(max(abs([kx_padded,ky_padded,kz_padded])));

for illumination = 1:n_illumination
    ki = k_incident(illumination,:);
    qx = KX-ki(1);
    qy = KY-ki(2);
    qz = KZ-ki(3);
    inside_fft = qx >= kx_padded(1)-interpolation_tolerance & ...
        qx <= kx_padded(end)+interpolation_tolerance & ...
        qy >= ky_padded(1)-interpolation_tolerance & ...
        qy <= ky_padded(end)+interpolation_tolerance & ...
        qz >= kz_padded(1)-interpolation_tolerance & ...
        qz <= kz_padded(end)+interpolation_tolerance;
    if any(propagating & ~inside_fft, 'all')
        error('born_rytov_voxel:InsufficientVoxelBandwidth', ...
            ['Decrease dx to resolve the full pre-objective Ewald hemisphere ', ...
             'for this incident direction.']);
    end

    query_x = min(max(qx(propagating),kx_padded(1)),kx_padded(end));
    query_y = min(max(qy(propagating),ky_padded(1)),ky_padded(end));
    query_z = min(max(qz(propagating),kz_padded(1)),kz_padded(end));
    sampled_object = complex(zeros(size(KX)));
    % Linear interpolation; add padding convergence or a NUFFT
    % when reconstruction accuracy, rather than portfolio scale, requires it.
    sampled_object(propagating) = interpn(kx_padded, ky_padded, kz_padded, F, ...
        query_x, query_y, query_z, 'linear');
    spectrum = propagator .* sampled_object;
    linear_field = fftshift(ifft2(ifftshift(spectrum))) / dx^2;
    detector_spectrum = spectrum .* pupil;
    born_spectrum(:,:,illumination) = detector_spectrum;
    uB(:,:,illumination) = ...
        fftshift(ifft2(ifftshift(detector_spectrum))) / dx^2;
    u0(:,:,illumination) = exp(1i * ...
        (ki(1)*X + ki(2)*Y + ki(3)*z_det));
    phiR(:,:,illumination) = ...
        linear_field ./ u0(:,:,illumination);
    rytov_field = u0(:,:,illumination) .* expm1(phiR(:,:,illumination));
    rytov_spectrum = fftshift(fft2(ifftshift(rytov_field))) * dx^2;
    rytov_spectrum(~pupil) = 0;
    uR(:,:,illumination) = ...
        fftshift(ifft2(ifftshift(rytov_spectrum))) / dx^2;
end

if any(~isfinite(uB), 'all') || any(~isfinite(uR), 'all')
    error('born_rytov_voxel:NonfiniteOutput', ...
        'The Born or Rytov field contains a non-finite value.');
end

info = struct('k0',k0, 'k',k, 'f_norm',norm(f(:)), ...
    'grid_size',grid_size, 'dx',dx, 'x',x, 'y',y, 'z',z, ...
    'kx',kx, 'ky',ky, 'kz',KZ, 'pupil',pupil, 'NA_det',NA_det, ...
    'k_incident',k_incident, 'z_det',z_det, ...
    'pad_factor',double(opts.pad_factor), 'interpolation','linear', ...
    'born_spectrum',born_spectrum, 'rytov_phase',phiR, ...
    'time',toc(started));
end
