function [uB, uR, u0, info] = born_rytov_spheroid( ...
    lambda, n_p, n_m, a, b, theta_inc, phi_inc, X, Y, Z, opts)
%BORN_RYTOV_SPHEROID Scalar Born and Rytov scattered fields of a spheroid.
% Angles are radians; a is the z semi-axis and b the equatorial radius.
% Both phase-free fields must converge independently within opts.tol.
% info.err is the larger normalized pair difference; failures use opts.on_fail.

started = tic;
if nargin < 11
    opts = struct;
end
positive_values = {lambda, n_m, a, b};
positive_names = {'lambda', 'n_m', 'a', 'b'};
for j = 1:numel(positive_values)
    validateattributes(positive_values{j}, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, positive_names{j});
end
validateattributes(n_p, {'numeric'}, {'finite', 'scalar'}, mfilename, 'n_p');
validateattributes(theta_inc, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'theta_inc');
validateattributes(phi_inc, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'phi_inc');
validateattributes(X, {'numeric'}, {'real', 'finite'}, mfilename, 'X');
validateattributes(Y, {'numeric'}, {'real', 'finite'}, mfilename, 'Y');
validateattributes(Z, {'numeric'}, {'real', 'finite'}, mfilename, 'Z');
if ~isequal(size(X), size(Y), size(Z))
    error('born_rytov_spheroid:InputSizeMismatch', 'X, Y, and Z must have equal sizes.');
end
if ~isstruct(opts) || ~isscalar(opts)
    error('born_rytov_spheroid:InvalidOptions', 'opts must be a scalar struct.');
end
known = {'tol', 'N_theta', 'N_max', 'on_fail', 'chunk'};
if any(~ismember(fieldnames(opts), known))
    error('born_rytov_spheroid:UnknownOption', 'opts contains an unknown field.');
end
lambda = double(lambda);
n_p = double(n_p);
n_m = double(n_m);
a = double(a);
b = double(b);
k = 2*pi*n_m/lambda;
f = k^2*((n_p/n_m)^2 - 1);
if ~isfinite(k) || k <= 0 || ~isfinite(f)
    error('born_rytov_spheroid:InvalidDerivedParameters', 'Derived k and f must be finite, with k positive.');
end
defaults = struct('tol', 1e-8, 'N_theta', ceil(2*k*max(a,b))+16, ...
    'N_max', 1200, 'on_fail', 'error', 'chunk', 2e6);
for j = 1:numel(known)
    name = known{j};
    if ~isfield(opts, name)
        opts.(name) = defaults.(name);
    end
end
validateattributes(opts.tol, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'opts.tol');
if opts.tol < 100*eps
    error('born_rytov_spheroid:InvalidTolerance', 'opts.tol must be at least 100*eps.');
end
for name = {'N_theta', 'N_max', 'chunk'}
    validateattributes(opts.(name{1}), {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', 'positive'}, mfilename, ['opts.' name{1}]);
end
if ~(ischar(opts.on_fail) && isrow(opts.on_fail) || ...
        isstring(opts.on_fail) && isscalar(opts.on_fail)) || ...
        ~any(strcmp(opts.on_fail, {'error', 'warn'}))
    error('born_rytov_spheroid:InvalidFailurePolicy', 'opts.on_fail must be error or warn.');
end
if ceil(1.5*double(opts.N_theta)) > double(opts.N_max)
    error('born_rytov_spheroid:InvalidQuadratureLimit', ...
        'The initial check order ceil(1.5*N_theta) must not exceed N_max.');
end
theta_inc = double(theta_inc);
phi_inc = double(phi_inc);
khat = [sin(theta_inc)*cos(phi_inc), sin(theta_inc)*sin(phi_inc), cos(theta_inc)];
shape = size(X);
pts = [double(X(:)), double(Y(:)), double(Z(:))];
rho_e = hypot(hypot(pts(:,1)/b, pts(:,2)/b), pts(:,3)/a);
if any(abs(rho_e-1) < 1e-9)
    error('born_rytov_spheroid:SurfacePoint', 'Please offset the point from the surface.');
end
inside = rho_e < 1;
rho = sqrt(sum(pts.^2, 2));
u0f = exp(1i*k*(pts*khat.'-rho));
V = 4*pi*a*b^2/3;
natural_scale = abs(f)*V./(4*pi*max(rho,max(a,b)));
phasefree = complex(zeros(size(inside)));
rytov_phasefree = complex(zeros(size(inside)));
orders = zeros(size(inside));
err = inf(size(inside));
unconverged = true(size(inside));
active = (1:size(pts,1)).';
N = double(opts.N_theta);
check_order = ceil(1.5*N);
% ponytail: global orders per failed group; use local subdivision for near-surface fields.
while ~isempty(active)
    pair_orders = [N, check_order];
    pair = complex(zeros(numel(active),2));
    for j = 1:2
        nodes = surface_nodes(a, b, k, khat, pair_orders(j));
        chunk_size = max(1, floor(double(opts.chunk)/size(nodes.xyz,1)));
        % ponytail: serial point chunks; add parfor if measured throughput needs it.
        for first = 1:chunk_size:numel(active)
            local = first:min(first+chunk_size-1,numel(active));
            index = active(local);
            pair(local,j) = born_eval(nodes,pts(index,:),f,k,khat,inside(index));
        end
    end
    scale = max(abs(pair(:,2)),natural_scale(active));
    difference = abs(pair(:,2)-pair(:,1));
    normalized = zeros(size(scale));
    normalized(scale>0) = difference(scale>0)./scale(scale>0);
    rytov_pair = u0f(active).*expm1(pair./u0f(active));
    rytov_scale = max(abs(rytov_pair(:,2)),natural_scale(active));
    rytov_difference = abs(rytov_pair(:,2)-rytov_pair(:,1));
    rytov_normalized = zeros(size(rytov_scale));
    rytov_normalized(rytov_scale>0) = ...
        rytov_difference(rytov_scale>0)./rytov_scale(rytov_scale>0);
    finite_pair = all(isfinite(pair),2) & all(isfinite(rytov_pair),2);
    phasefree(active) = pair(:,2);
    rytov_phasefree(active) = rytov_pair(:,2);
    orders(active) = pair_orders(2);
    normalized = max(normalized,rytov_normalized);
    normalized(~finite_pair) = inf;
    err(active) = normalized;
    unconverged(active) = ~finite_pair | ~(difference <= double(opts.tol)*scale & ...
        rytov_difference <= double(opts.tol)*rytov_scale);
    active = active(unconverged(active));
    if check_order == double(opts.N_max)
        break
    end
    if ceil(1.5*(2*N)) <= double(opts.N_max)
        N = 2*N;
        check_order = ceil(1.5*N);
    else
        cap_base = floor(double(opts.N_max)/1.5);
        N = max(N,cap_base);
        check_order = double(opts.N_max);
    end
end
% ponytail: double common-phase reduction; higher precision beyond the tested range.
common_phase = exp(1i*k*rho);
u0 = reshape(u0f .* common_phase, shape);
uB = reshape(phasefree .* common_phase, shape);
uR = reshape(rytov_phasefree .* common_phase, shape);
nonfinite = ~isfinite(uB(:)) | ~isfinite(uR(:));
unconverged = unconverged | nonfinite;
err(nonfinite) = inf;
info = struct('k', k, 'f', f, 'khat', khat, 'V', V, ...
    'N_theta', reshape(orders,shape), 'err', reshape(err,shape), ...
    'inside', reshape(inside, shape), 'unconverged', reshape(unconverged,shape), ...
    'validated', ~any(unconverged), 'time', toc(started));
if any(unconverged)
    message = sprintf('%d point(s) failed quadrature convergence or produced non-finite output; minimum abs(rho_e-1) = %.6g.', ...
        nnz(unconverged),min(abs(rho_e(unconverged)-1)));
    if strcmp(opts.on_fail,'error')
        error('born_rytov_spheroid:ConvergenceFailure','%s',message);
    end
    warning('born_rytov_spheroid:ConvergenceFailure','%s',message);
end
end

function nodes = surface_nodes(a, b, k, khat, N_theta)
[g, weights] = gauss_legendre(N_theta);
theta = pi*(g+1)/2;
phi = (0:2*N_theta-1) * (pi/N_theta);
[T, P] = ndgrid(theta, phi);
area_weight = (pi*weights/2) * ones(1, 2*N_theta) * (pi/N_theta);
s = sin(T(:));
c = cos(T(:));
nodes.xyz = [b*s.*cos(P(:)), b*s.*sin(P(:)), a*c];
nodes.normal = [a*s.*cos(P(:)), a*s.*sin(P(:)), b*c] .* (b*s.*area_weight(:));
projection = nodes.xyz*khat.';
incident = exp(1i*k*projection);
nodes.w = projection .* incident / (2i*k);
nodes.dnw = (nodes.normal*khat.') .* incident .* (1+1i*k*projection) / (2i*k);
nodes.radius2 = sum(nodes.xyz.^2, 2);
end

function value = born_eval(nodes, pts, f, k, khat, inside)
% Return the integral with its common exp(ik*|r|) phase removed.
dx = nodes.xyz(:,1).' - pts(:,1);
dy = nodes.xyz(:,2).' - pts(:,2);
dz = nodes.xyz(:,3).' - pts(:,3);
R = sqrt(dx.^2 + dy.^2 + dz.^2);
rho = sqrt(sum(pts.^2, 2));
deltaR = (nodes.radius2.' - 2*(pts*nodes.xyz.')) ./ (R+rho);
G = exp(1i*k*deltaR) ./ (4*pi*R);
normal_distance = dx.*nodes.normal(:,1).' + ...
    dy.*nodes.normal(:,2).' + dz.*nodes.normal(:,3).';
dnG = (1i*k - 1./R) .* G .* normal_distance ./ R;
value = G*nodes.dnw - dnG*nodes.w;
projection = pts*khat.';
jump = projection .* exp(1i*k*(projection-rho)) / (2i*k);
value = f * (value - inside.*jump);
end
