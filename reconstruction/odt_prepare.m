function data = odt_prepare(U, u0, geom, model, opts)
%ODT_PREPARE Physical detector fields -> nonuniform Fourier samples of chi.
% chi = n^2/n_m^2 - 1 is REAL (nonabsorbing, isotropic sample).
% geom is the info returned by born-rytov/born_rytov_voxel. U and u0 have
% size NdetX-by-NdetY-by-Nillum. model is 'born' or 'rytov'.
% geom.grid_size and geom.dx specify the reconstruction volume; optional
% geom.detector_dx specifies a different detector pitch (default geom.dx).
% The measured Rytov logarithm is taken AFTER the physical detector pupil;
% its linear Fourier interpretation remains an approximation.
% opts.unwrapped_phase optionally supplies angle(U/u0) on its correct branch.
% opts.phase_floor (1e-6) bounds abs(U/u0) away from the Rytov log singularity.
% No ground-truth volume, spectrum, or phase in geom is read.

if nargin < 5, opts = struct; end
model = validatestring(model, {'born','rytov'});
validateattributes(U, {'numeric'}, {'nonempty','finite'}, mfilename,'U');
validateattributes(u0, {'numeric'}, {'nonempty','finite'}, mfilename,'u0');
if ~isstruct(geom) || ~isscalar(geom) || ...
        ~all(isfield(geom, {'grid_size','dx','k','k0','kx','ky','kz', ...
        'pupil','k_incident','z_det','NA_det'}))
    error('odt_prepare:InvalidGeometry', 'Use geometry from born_rytov_voxel.');
end
if ~isstruct(opts) || ~isscalar(opts) || ...
        any(~ismember(fieldnames(opts),{'unwrapped_phase','phase_floor','padding','weights'}))
    error('odt_prepare:InvalidOptions', 'Unknown or invalid preparation options.');
end
if ~isfield(opts,'phase_floor'), opts.phase_floor = 1e-6; end
validateattributes(opts.phase_floor, {'numeric'}, ...
    {'real','finite','scalar','positive','<',1});
sz = double(geom.grid_size(:).');
validateattributes(sz, {'numeric'}, {'integer','finite','numel',3,'>=',4});
validateattributes(geom.dx, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(geom.k, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(geom.k0, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(geom.NA_det, {'numeric'}, ...
    {'real','finite','scalar','positive','<',geom.k/geom.k0});
validateattributes(geom.z_det, {'numeric'}, ...
    {'real','finite','scalar','>',(sz(3)/2-0.5)*geom.dx});
validateattributes(geom.k_incident, {'numeric'}, ...
    {'real','finite','ncols',3,'nonempty'});
ni = size(geom.k_incident,1);
det_size=[numel(geom.kx),numel(geom.ky)];
det_dx=geom.dx;
if isfield(geom,'detector_dx'), det_dx=geom.detector_dx; end
validateattributes(det_dx,{'numeric'},{'real','finite','scalar','positive'});
if any(mod(sz,2)) || ~isequal(size(U),size(u0)) || ...
        any(mod(det_size,2)) || ...
        ~isequal([size(U,1),size(U,2),size(U,3)],[det_size,ni]) || ndims(U)>3
    error('odt_prepare:InvalidSize', 'Field and geometry dimensions must agree.');
end
if any(abs(u0(:)) <= realmin)
    error('odt_prepare:InvalidReference', 'Reference field must be nonzero.');
end
dx = geom.dx; k = geom.k;
[KX,KY] = ndgrid(geom.kx,geom.ky);
pupil = KX.^2+KY.^2 <= (geom.k0*geom.NA_det)^2;
if ~isequal(size(pupil),det_size) || ~isequal(pupil,geom.pupil) || ~any(pupil(:))
    error('odt_prepare:InvalidGeometry', 'Detector pupil does not match geometry.');
end
kz = sqrt(k^2-KX(pupil).^2-KY(pupil).^2);
ki = geom.k_incident;
if any(ki(:,3)<=0) || any(abs(sum(ki.^2,2)-k^2)>1e-10*k^2) || ...
        any(hypot(ki(:,1),ki(:,2))>geom.k0*geom.NA_det+64*eps(k)) || ...
        any(abs(kz-geom.kz(pupil))>1e-10*k)
    error('odt_prepare:InvalidGeometry', 'Illumination or detector kz is inconsistent.');
end
if strcmp(model,'born') && isfield(opts,'unwrapped_phase')
    error('odt_prepare:InvalidOptions', 'unwrapped_phase applies only to Rytov.');
end
if isfield(opts,'unwrapped_phase')
    validateattributes(opts.unwrapped_phase, {'numeric'}, {'real','finite'});
    if ~isequal(size(opts.unwrapped_phase),size(U))
        error('odt_prepare:InvalidPhaseSize', 'Unwrapped phase must have the size of U.');
    end
end
[X,Y] = ndgrid((-det_size(1)/2:det_size(1)/2-1)*det_dx, ...
    (-det_size(2)/2:det_size(2)/2-1)*det_dx);
ns = nnz(pupil);
q = zeros(ns*ni,3); samples = complex(zeros(ns*ni,1));
phase_order_error = zeros(ni,1);
border_phase = zeros(ni,1);
for j = 1:ni
    ref = exp(1i*(ki(j,1)*X+ki(j,2)*Y+ki(j,3)*geom.z_det));
    actual_ref = u0(:,:,j);
    amplitude = mean(actual_ref(:)./ref(:));
    if abs(amplitude)<=realmin || ...
            norm(actual_ref(:)-amplitude*ref(:))>1e-8*norm(actual_ref(:))
        error('odt_prepare:InvalidReference', ...
            'Reference must be the stated plane wave, up to one complex amplitude per angle.');
    end
    if strcmp(model,'born')
        linear_field = (U(:,:,j)-actual_ref)/amplitude;
    else
        ratio = U(:,:,j)./actual_ref;
        if any(~isfinite(ratio(:))) || any(abs(ratio(:))<=opts.phase_floor)
            error('odt_prepare:PhaseSingularity', 'Rytov log is undefined or poorly conditioned.');
        end
        wrapped = angle(ratio);
        if isfield(opts,'unwrapped_phase')
            phase = opts.unwrapped_phase(:,:,j);
            if max(abs(angle(exp(1i*(phase(:)-wrapped(:)))))) > 1e-6
                error('odt_prepare:InconsistentPhase', 'Supplied phase does not match U/u0 modulo 2*pi.');
            end
        else
            [phase,phase_order_error(j),border_phase(j)] = unwrap_checked(wrapped);
        end
        linear_field = ref.*(log(abs(ratio))+1i*phase);
    end
    field_ft = fftshift(fft2(ifftshift(linear_field)))*det_dx^2;
    rows = (j-1)*ns+(1:ns);
    % Invert i*exp(i*kz*z_det)/(2*kz); f=k^2*chi in the sibling forward.
    samples(rows) = field_ft(pupil).*(-2i*kz).*exp(-1i*kz*geom.z_det)/k^2;
    q(rows,:) = [KX(pupil)-ki(j,1),KY(pupil)-ki(j,2),kz-ki(j,3)];
end
if any(~isfinite(samples))
    error('odt_prepare:NonfiniteData', 'Transformed Fourier data are non-finite.');
end
if any(abs(q)>=pi/dx,'all')
    error('odt_prepare:InsufficientBandwidth', 'Decrease dx to represent all Ewald samples.');
end
padding=2;
if isfield(geom,'pad_factor'), padding=geom.pad_factor; end
if isfield(opts,'padding'), padding=opts.padding; end
weights=ones(ni,1);
if isfield(opts,'weights'), weights=opts.weights; end
validateattributes(weights,{'numeric'},{'real','finite','positive','vector','numel',ni});
op=odt_operator(q,sz,dx,padding,repelem(weights(:),ns));
y=samples.*op.sample_scale;
data = struct('operator',op,'y',y, ...
    'n_m',k/geom.k0,'dx',dx,'model',model,'NA_det',geom.NA_det, ...
    'grid_size',sz,'detector_size',det_size,'detector_dx',det_dx, ...
    'q',q,'samples',samples,'weights',weights(:),'padding',padding, ...
    'phase_order_error',phase_order_error,'border_phase',border_phase, ...
    'scalar_validity','requires independent validation','gridding','matched trilinear', ...
    'symmetry','Hermitian symmetry inferred from real chi');
end

function [phase,disagreement,border_value] = unwrap_checked(wrapped)
% Two path orders are a consistency test, not a sampling or branch proof.
border = false(size(wrapped)); border([1,end],:) = true; border(:,[1,end]) = true;
phase = unwrap(unwrap(wrapped,[],1),[],2);
other = unwrap(unwrap(wrapped,[],2),[],1);
phase = phase-2*pi*round(median(phase(border))/(2*pi));
other = other-2*pi*round(median(other(border))/(2*pi));
disagreement = max(abs(phase(:)-other(:)));
border_value = max(abs(phase(border)));
if disagreement>1e-6 || border_value>=pi/2
    error('odt_prepare:AmbiguousPhase', ...
        ['Automatic phase unwrapping failed its path/border checks. ', ...
        'Provide opts.unwrapped_phase or increase spatial sampling/FOV.']);
end
end
