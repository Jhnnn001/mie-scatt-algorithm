function [Ex, Ey, Ez, info] = spheroid_field(lambda, n_p, n_m, a, b, ...
    theta_inc, phi_inc, pol, X, Y, Z, opts)
%SPHEROID_FIELD Electric field of a homogeneous spheroid.

if nargin < 12 || isempty(opts)
    opts = struct();
end
[opts, field, side] = checked_inputs(lambda, n_p, n_m, a, b, ...
    theta_inc, phi_inc, pol, X, Y, Z, opts);
started = tic;
k_m = 2*pi*n_m/lambda;
k_p = 2*pi*n_p/lambda;
inside_envelope = envelope_member(k_m, k_p, n_p/n_m, a, b);

if abs(a - b) <= 64*eps(max(a, b))
    warn_if_outside(inside_envelope);
    [E, ~] = mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
        pol, X, Y, Z, field, side);
    mie_tolerance = struct('computed', true, 'value', 1e-10, ...
        'limit', opts.tol, 'passed', opts.tol >= 1e-10);
    info = route_info('sphere', 'mie', Inf, 0, 0, ...
        inside_envelope && mie_tolerance.passed, ...
        struct('mie_solution', true, 'mie_tolerance', mie_tolerance), ...
        opts, started);
    info.sol = struct('route', 'mie', 'radius', a);
    [Ex, Ey, Ez] = split_field(E, size(X));
    return
end

if a > b
    type = 'prolate';
else
    type = 'oblate';
end
semifocal = sph_semifocal(a, b);
xi0 = a/semifocal;

if n_p == n_m
    warn_if_outside(inside_envelope);
    enforce_side(type, semifocal, xi0, X, Y, Z, side);
    [E, ~] = incident_plane_wave(k_m, n_m, theta_inc, phi_inc, ...
        pol, X, Y, Z);
    if strcmp(field, 'scattered')
        E(:) = 0;
    end
    checks = struct('analytic_identity', true);
    info = route_info(type, 'zero_contrast', xi0, ...
        k_m*semifocal, k_p*semifocal, inside_envelope, checks, ...
        opts, started);
    info.sol = struct('route', 'zero_contrast', 'type', type, ...
        'semifocal', semifocal, 'xi0', xi0);
    [Ex, Ey, Ez] = split_field(E, size(X));
    return
end

if strcmp(field, 'incident')
    warn_if_outside(inside_envelope);
    enforce_side(type, semifocal, xi0, X, Y, Z, side);
    [E, ~] = incident_plane_wave(k_m, n_m, theta_inc, phi_inc, ...
        pol, X, Y, Z);
    checks = struct('analytic_incident', true);
    info = route_info(type, 'incident', xi0, ...
        k_m*semifocal, k_p*semifocal, inside_envelope, checks, ...
        opts, started);
    info.sol = struct('route', 'incident', 'type', type, ...
        'semifocal', semifocal, 'xi0', xi0);
    [Ex, Ey, Ez] = split_field(E, size(X));
    return
end

sol = spheroid_solve(lambda, n_p, n_m, a, b, theta_inc, ...
    phi_inc, pol, opts);
[E, ~] = spheroid_eval(sol, X, Y, Z, field, side);
if isfield(sol, 'info')
    info = sol.info;
else
    info = struct();
end
info.sol = sol;
info.type = type;
info.xi0 = xi0;
info.c_ext = k_m*semifocal;
info.c_int = k_p*semifocal;
info.time = toc(started);
[Ex, Ey, Ez] = split_field(E, size(X));
end

function warn_if_outside(inside_envelope)
if ~inside_envelope
    warning('spheroid_field:OutsideValidatedEnvelope', ...
        'The solution is outside the tested parameter envelope.');
end
end

function [opts, field, side] = checked_inputs(lambda, n_p, n_m, a, b, ...
    theta_inc, phi_inc, pol, X, Y, Z, opts)
validateattributes(lambda, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'lambda');
validateattributes(n_p, {'numeric'}, {'finite', 'scalar'}, mfilename, 'n_p');
validateattributes(n_m, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'n_m');
validateattributes(a, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'a');
validateattributes(b, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'b');
validateattributes(theta_inc, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'theta_inc');
validateattributes(phi_inc, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'phi_inc');
validateattributes(pol, {'numeric'}, ...
    {'finite', 'vector', 'numel', 2}, mfilename, 'pol');
validateattributes(X, {'numeric'}, {'real', 'finite'}, mfilename, 'X');
validateattributes(Y, {'numeric'}, {'real', 'finite'}, mfilename, 'Y');
validateattributes(Z, {'numeric'}, {'real', 'finite'}, mfilename, 'Z');
if real(n_p) <= 0 || imag(n_p) < 0
    error('spheroid_field:InvalidRefractiveIndex', ...
        'n_p must have positive real part and nonnegative imaginary part.');
end
if ~any(pol ~= 0)
    error('spheroid_field:InvalidPolarization', ...
        'pol must contain a nonzero component.');
end
if ~isequal(size(X), size(Y), size(Z))
    error('spheroid_field:InputSizeMismatch', ...
        'X, Y, and Z must have equal sizes.');
end
if ~isstruct(opts) || ~isscalar(opts)
    error('spheroid_field:InvalidOptions', 'opts must be a scalar structure.');
end
defaults = struct('tol', 1e-6, 'field', 'total', 'side', 'auto', ...
    'Lmax', 400, 'Mmax', 300, 'on_fail', 'error', 'verbose', false);
names = fieldnames(opts);
allowed = fieldnames(defaults);
if any(~ismember(names, allowed))
    error('spheroid_field:InvalidOptions', ...
        'opts contains an unknown field.');
end
for k = 1:numel(allowed)
    name = allowed{k};
    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end
validateattributes(opts.tol, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'opts.tol');
validateattributes(opts.Lmax, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', '>=', 1, '<=', 400}, ...
    mfilename, 'opts.Lmax');
validateattributes(opts.Mmax, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', '>=', 1, '<=', 300}, ...
    mfilename, 'opts.Mmax');
if ~(islogical(opts.verbose) && isscalar(opts.verbose))
    error('spheroid_field:InvalidOptions', ...
        'opts.verbose must be a logical scalar.');
end
field = checked_option(opts.field, ...
    {'total', 'scattered', 'incident'});
side = checked_option(opts.side, {'auto', 'in', 'out'});
opts.on_fail = checked_option(opts.on_fail, {'error', 'warn'});
opts.field = field;
opts.side = side;
end

function value = checked_option(value, choices)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~(ischar(value) && isrow(value) && ismember(value, choices))
    error('spheroid_field:InvalidOption', ...
        'Option must be one of: %s.', strjoin(choices, ', '));
end
end

function enforce_side(type, semifocal, xi0, X, Y, Z, side)
if strcmp(side, 'auto')
    return
end
g = sph_coords(type, X(:)/semifocal, Y(:)/semifocal, Z(:)/semifocal);
if any(abs(g.xi - xi0) > 1e-8)
    error('spheroid_field:InvalidForcedSide', ...
        'Forced side is allowed only within 1e-8 in xi of the boundary.');
end
end

function yes = envelope_member(k_m, k_p, n_relative, a, b)
scale = max(a, b);
aspect = a/b;
loss_ratio = imag(n_relative)/real(n_relative);
x_ext = k_m*scale;
yes = x_ext >= 0.02 && x_ext <= 250 && abs(k_p)*scale <= 250 && ...
    aspect >= 1/4 && aspect <= 4 && real(n_relative) >= 0.5 && ...
    real(n_relative) <= 2 && loss_ratio <= 0.1;
end

function info = route_info(type, route, xi0, c_ext, c_int, validated, ...
    checks, opts, started)
info = struct('type', type, 'route', route, 'xi0', xi0, ...
    'c_ext', c_ext, 'c_int', c_int, 'checks', checks, ...
    'validated', validated, 'options', opts, 'time', toc(started));
end

function [Ex, Ey, Ez] = split_field(E, shape)
Ex = reshape(E(:, 1), shape);
Ey = reshape(E(:, 2), shape);
Ez = reshape(E(:, 3), shape);
end
