function sol = spheroid_solve(lambda, n_p, n_m, a, b, theta_inc, ...
    phi_inc, pol, opts)
%SPHEROID_SOLVE Solve plane-wave scattering by a homogeneous spheroid.

narginchk(8, 9);
if nargin < 9 || isempty(opts)
    opts = struct();
end
opts = checked_inputs(lambda, n_p, n_m, a, b, theta_inc, ...
    phi_inc, pol, opts);
if abs(a - b) <= 64*eps(max(a, b))
    error('spheroid_solve:SphereRouteRequired', ...
        'Use spheroid_field for geometries on the Mie sphere route.');
end
if n_p == n_m
    error('spheroid_solve:ZeroContrastRouteRequired', ...
        'Use spheroid_field for the analytic zero-contrast route.');
end

problem = make_problem(lambda, n_p, n_m, a, b, theta_inc, phi_inc, pol);
problem.axial_incidence = abs(sin(theta_inc)) <= min(8*eps, opts.tol);
problem.opts = opts;
[L, M, Qtheta, Qphi] = starting_sizes(problem, opts);
points = deterministic_points(problem);
checks = empty_checks();
last_solution = [];
failure_reason = '';
completed = false;

for round_index = 1:5
    try
        base = solve_once(problem, L, M, Qtheta, Qphi, {});
    catch exception
        if isempty(last_solution)
            rethrow(exception)
        end
        failure_reason = sprintf('base solve failed: %s', exception.message);
        break
    end
    last_solution = base;
    checks = empty_checks();
    checks.special_functions = special_function_check(base.mu_data);
    checks.incident_tail = incident_tail_check(base);
    try
        checks.qphi_independence = qphi_check(base, 2*Qphi);
    catch exception
        checks.qphi_independence = failed_check(exception);
    end
    try
        checks.modal_tail = modal_tail_check(base);
    catch exception
        checks.modal_tail = failed_check(exception);
    end

    if ~checks.special_functions.passed
        failure_reason = 'a special-function verification failed';
        break
    end
    if ~checks.incident_tail.passed
        [M, L, Qtheta, Qphi, possible] = grow_M( ...
            M, L, Qtheta, Qphi, problem, opts);
        if ~possible
            failure_reason = 'the incident spectrum reached Mmax';
            break
        end
        continue
    end
    if ~checks.qphi_independence.passed
        if ~isfinite(checks.qphi_independence.value)
            failure_reason = 'the Qphi independence check failed';
            break
        end
        Qphi = 2*Qphi;
        continue
    end
    if ~checks.modal_tail.passed
        if ~isfinite(checks.modal_tail.value)
            failure_reason = 'the modal-tail check failed';
            break
        end
        [L, Qtheta, possible] = grow_L(L, Qtheta, problem, opts);
        if ~possible
            failure_reason = 'the modal tail reached Lmax';
            break
        end
        continue
    end

    try
        checks.surface_residual = surface_residual_check(base);
    catch exception
        checks.surface_residual = failed_check(exception);
    end
    if ~checks.surface_residual.passed
        if ~isfinite(checks.surface_residual.value)
            failure_reason = 'the independent surface check failed';
            break
        end
        [L, Qtheta, possible_L] = grow_L(L, Qtheta, problem, opts);
        [M, L, Qtheta, Qphi, possible_M] = grow_M( ...
            M, L, Qtheta, Qphi, problem, opts);
        Qtheta = ceil(1.5*Qtheta);
        if ~(possible_L && possible_M)
            failure_reason = 'the surface check reached a cutoff ceiling';
            break
        end
        continue
    end

    refined_Qtheta = ceil(1.5*Qtheta);
    try
        base_signature = solution_signature(base, points);
        candidate = solve_once(problem, L, M, refined_Qtheta, Qphi, ...
            base.mu_data);
        checks.qtheta_independence = comparison_check(base, candidate, ...
            points, base_signature);
    catch exception
        checks.qtheta_independence = failed_check(exception);
    end
    if ~checks.qtheta_independence.passed
        if ~isfinite(checks.qtheta_independence.value)
            failure_reason = 'the Qtheta verification solve failed';
            break
        end
        Qtheta = refined_Qtheta;
        continue
    end

    cutoff_L = L + 5;
    cutoff_M = M + 2;
    if cutoff_L > opts.Lmax || cutoff_M > opts.Mmax || cutoff_M > cutoff_L
        checks.cutoff_stability = make_check(false, Inf, false, ...
            struct('reason', 'cutoff ceiling'));
        failure_reason = 'the cutoff verification exceeds Lmax or Mmax';
        break
    end
    cutoff_Qtheta = max(Qtheta, minimum_Qtheta(problem, cutoff_L));
    cutoff_Qphi = max(Qphi, minimum_Qphi(problem, cutoff_M));
    try
        candidate = solve_once(problem, cutoff_L, cutoff_M, ...
            cutoff_Qtheta, cutoff_Qphi, {});
        checks.cutoff_stability = comparison_check(base, candidate, ...
            points, base_signature);
    catch exception
        checks.cutoff_stability = failed_check(exception);
    end
    if ~checks.cutoff_stability.passed
        if ~isfinite(checks.cutoff_stability.value)
            failure_reason = 'the cutoff verification solve failed';
            break
        end
        L = cutoff_L;
        M = cutoff_M;
        Qtheta = cutoff_Qtheta;
        Qphi = cutoff_Qphi;
        continue
    end

    completed = all_checks_passed(checks);
    if completed
        break
    end
end

if isempty(last_solution)
    error('spheroid_solve:NoSolution', ...
        'No finite spheroidal solution was produced.');
end
if ~completed && isempty(failure_reason)
    failure_reason = 'the five-round convergence limit was reached';
end
inside_envelope = envelope_member(problem);
last_solution.info = make_info(last_solution, checks, round_index, ...
    inside_envelope, completed, failure_reason, opts);
last_solution.opts = opts;
sol = last_solution;

if ~inside_envelope
    warning('spheroid_solve:OutsideValidatedEnvelope', ...
        'The solution is outside the tested parameter envelope.');
end
if ~completed
    message = sprintf('Spheroidal convergence was not verified: %s.', ...
        failure_reason);
    if strcmp(opts.on_fail, 'error')
        error('spheroid_solve:NoConvergence', '%s', message);
    end
    warning('spheroid_solve:NoConvergence', '%s', message);
end
end

function opts = checked_inputs(lambda, n_p, n_m, a, b, theta_inc, ...
    phi_inc, pol, opts)
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
if real(n_p) <= 0 || imag(n_p) < 0
    error('spheroid_solve:InvalidRefractiveIndex', ...
        'n_p must have positive real part and nonnegative imaginary part.');
end
if ~any(pol ~= 0)
    error('spheroid_solve:InvalidPolarization', ...
        'pol must contain a nonzero component.');
end
if ~isstruct(opts) || ~isscalar(opts)
    error('spheroid_solve:InvalidOptions', 'opts must be a scalar structure.');
end
defaults = struct('tol', 1e-6, 'field', 'total', 'side', 'auto', ...
    'Lmax', 400, 'Mmax', 300, 'on_fail', 'error', 'verbose', false);
names = fieldnames(opts);
allowed = fieldnames(defaults);
if any(~ismember(names, allowed))
    error('spheroid_solve:InvalidOptions', ...
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
    error('spheroid_solve:InvalidOptions', ...
        'opts.verbose must be a logical scalar.');
end
opts.field = checked_choice(opts.field, ...
    {'total', 'scattered', 'incident'});
opts.side = checked_choice(opts.side, {'auto', 'in', 'out'});
opts.on_fail = checked_choice(opts.on_fail, {'error', 'warn'});
end

function value = checked_choice(value, choices)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~(ischar(value) && isrow(value) && ismember(value, choices))
    error('spheroid_solve:InvalidOption', ...
        'Option must be one of: %s.', strjoin(choices, ', '));
end
end

function problem = make_problem(lambda, n_p, n_m, a, b, theta_inc, ...
    phi_inc, pol)
if a > b
    type = 'prolate';
else
    type = 'oblate';
end
semifocal = sph_semifocal(a, b);
k0 = 2*pi/lambda;
km = n_m*k0;
kp = n_p*k0;
[~, ~, ~, E0] = incident_plane_wave(km, n_m, theta_inc, ...
    phi_inc, pol, 0, 0, 0);
problem = struct('lambda', lambda, 'n_p', n_p, 'n_m', n_m, ...
    'a', a, 'b', b, 'theta_inc', theta_inc, 'phi_inc', phi_inc, ...
    'pol', pol(:).', 'k0', k0, 'km', km, 'kp', kp, ...
    'semifocal', semifocal, 'type', type, 'xi0', a/semifocal, ...
    'c_ext', km*semifocal, 'c_int', kp*semifocal, 'E0', E0, ...
    'E0_norm', norm(E0), 'chi', abs(km*b*sin(theta_inc)));
end

function [L, M, Qtheta, Qphi] = starting_sizes(problem, opts)
x = max(problem.km*max(problem.a, problem.b), ...
    abs(problem.kp)*max(problem.a, problem.b));
L = ceil(x + 4.05*x^(1/3) + 10);
M = ceil(problem.chi + 4*problem.chi^(1/3) + 10) + 1;
M = max(1, min(M, L));
if L > opts.Lmax || M > opts.Mmax
    error('spheroid_solve:CutoffLimit', ...
        'Starting cutoffs L=%d, M=%d exceed Lmax or Mmax.', L, M);
end
Qtheta = minimum_Qtheta(problem, L);
Qphi = minimum_Qphi(problem, M);
end

function value = minimum_Qtheta(problem, L)
value = ceil(L + max(problem.c_ext, abs(problem.c_int)) + 30);
end

function value = minimum_Qphi(problem, M)
value = ceil(2*(M + problem.chi) + 32);
end

function sol = solve_once(problem, L, M, Qtheta, Qphi, cache)
if M > L || L > 400 || M > 300
    error('spheroid_solve:CutoffLimit', 'The requested cutoffs are invalid.');
end
[g_surface, theta, weights] = surface_meridian(problem, Qtheta);
incident = incident_spectrum(problem, theta, Qphi);
if isempty(cache)
    mu_data = cell(M + 1, 1);
else
    valid_cache = numel(cache) == M + 1;
    if valid_cache
        for mu = 0:M
            item = cache{mu + 1};
            unused = problem.axial_incidence && mu ~= 1;
            valid_cache = valid_cache && ((unused && isempty(item)) || ...
                (~isempty(item) && item.eg_ext.L == L && item.eg_int.L == L));
        end
    end
    if ~valid_cache
        error('spheroid_solve:InvalidCache', ...
            'A solve cache must match both L and M.');
    end
    mu_data = cache;
end
mode_template = struct('m', 0, 'mu', 0, 'a', [], 'b', [], ...
    'c', [], 'd', [], 'rhs', [], 'rhs_norm', NaN, 'rcond', NaN, ...
    'numerical_rank', NaN, 'solver_method', '', ...
    'analytically_zeroed', false, 'roundoff_zeroed', false);
modes = repmat(mode_template, 2*M + 1, 1);
for mu = 0:M
    count = L - mu + 1;
    if problem.axial_incidence && mu ~= 1
        for m = signed_modes(mu)
            k = m + M + 1;
            zero = complex(zeros(count, 1));
            modes(k) = struct('m', m, 'mu', mu, ...
                'a', zero, 'b', zero, 'c', zero, 'd', zero, ...
                'rhs', complex(zeros(4*count, 1)), 'rhs_norm', 0, ...
                'rcond', NaN, 'numerical_rank', NaN, ...
                'solver_method', 'analytic_zero', ...
                'analytically_zeroed', true, 'roundoff_zeroed', true);
        end
        continue
    end
    if isempty(mu_data{mu + 1})
        mu_data{mu + 1} = make_mu_data(problem, mu, L);
    end
    data = mu_data{mu + 1};
    Wext = sph_angular(data.eg_ext, g_surface);
    Wint = sph_angular(data.eg_int, g_surface);
    Vext = sph_radial_ode_eval(data.rod_ext, g_surface);
    [Vint, internal_ok] = sph_radial(data.eg_int, 1, ...
        g_surface, problem.c_int);
    if ~internal_ok
        error('spheroid_solve:RadialFailure', ...
            'The internal radial series failed for mu=%d.', mu);
    end
    for m = signed_modes(mu)
        [Mext, Next] = sph_vecwave(Wext, Vext, g_surface, m, ...
            problem.c_ext);
        [Mint, Nint] = sph_vecwave(Wint, Vint, g_surface, m, ...
            problem.c_int);
        [matrix, rhs] = assemble_mode(Wext.W0, weights, incident, ...
            Mext, Next, Mint, Nint, m, problem.n_m, problem.n_p);
        scale = max(abs(matrix), [], 1);
        if any(~isfinite(scale)) || any(scale <= 0)
            error('spheroid_solve:InvalidColumnScale', ...
                'A Galerkin column has zero or nonfinite scale for m=%d.', m);
        end
        equilibrated = matrix ./ scale;
        condition_reciprocal = rcond(equilibrated);
        rhs_norm = norm(rhs);
        axial_incidence = problem.axial_incidence;
        analytically_zeroed = all(rhs == 0) || ...
            (axial_incidence && abs(m) ~= 1);
        if analytically_zeroed
            coefficients = complex(zeros(size(rhs)));
            numerical_rank = NaN;
            solver_method = 'analytic_zero';
        elseif condition_reciprocal <= size(equilibrated, 1)*eps
            coefficients = lsqminnorm(equilibrated, rhs) ./ scale.';
            numerical_rank = rank(equilibrated);
            solver_method = 'lsqminnorm';
        else
            coefficients = (equilibrated \ rhs) ./ scale.';
            numerical_rank = size(equilibrated, 2);
            solver_method = 'backslash';
        end
        if any(~isfinite(real(coefficients))) || ...
                any(~isfinite(imag(coefficients)))
            error('spheroid_solve:NonfiniteCoefficients', ...
                'The Galerkin solve produced nonfinite coefficients.');
        end
        k = m + M + 1;
        modes(k) = struct('m', m, 'mu', mu, ...
            'a', coefficients(1:count), ...
            'b', coefficients(count + (1:count)), ...
            'c', coefficients(2*count + (1:count)), ...
            'd', coefficients(3*count + (1:count)), ...
            'rhs', rhs, 'rhs_norm', rhs_norm, ...
            'rcond', condition_reciprocal, ...
            'numerical_rank', numerical_rank, ...
            'solver_method', solver_method, ...
            'analytically_zeroed', analytically_zeroed, ...
            'roundoff_zeroed', analytically_zeroed);
    end
end
sol = problem;
sol.route = 'spheroid';
sol.L = L;
sol.M = M;
sol.Qtheta = Qtheta;
sol.Qphi = Qphi;
sol.mu_data = mu_data;
sol.modes = modes;
end

function values = signed_modes(mu)
if mu == 0
    values = 0;
else
    values = [-mu, mu];
end
end

function data = make_mu_data(problem, mu, L)
eg_ext = sph_eigen(problem.type, mu, problem.c_ext, L);
eg_int = sph_eigen(problem.type, mu, problem.c_int, L);
[rod_ext, radial_check] = sph_radial_ode(eg_ext, problem.c_ext, ...
    problem.xi0, 1.5);
if ~radial_check.ok
    error('spheroid_solve:RadialFailure', ...
        'Outgoing radial verification failed for mu=%d.', mu);
end
data = struct('mu', mu, 'eg_ext', eg_ext, 'eg_int', eg_int, ...
    'rod_ext', rod_ext, 'radial_check', radial_check);
end

function [matrix, rhs] = assemble_mode(Sext, weights, incident, ...
    Mext, Next, Mint, Nint, m, n_m, n_p)
project = @(F) Sext * (weights.' .* F).';
matrix = [project(Next.eta), project(Mext.eta), ...
          -project(Nint.eta), -project(Mint.eta); ...
          project(Next.phi), project(Mext.phi), ...
          -project(Nint.phi), -project(Mint.phi); ...
          -1i*n_m*project(Mext.eta), -1i*n_m*project(Next.eta), ...
          1i*n_p*project(Mint.eta), 1i*n_p*project(Nint.eta); ...
          -1i*n_m*project(Mext.phi), -1i*n_m*project(Next.phi), ...
          1i*n_p*project(Mint.phi), 1i*n_p*project(Nint.phi)];
bin = mod(m, incident.Qphi) + 1;
rhs = -[Sext*(weights.*incident.Eeta(:, bin)); ...
        Sext*(weights.*incident.Ephi(:, bin)); ...
        Sext*(weights.*incident.Heta(:, bin)); ...
        Sext*(weights.*incident.Hphi(:, bin))];
end

function incident = incident_spectrum(problem, theta, Qphi)
phi = 2*pi*(0:Qphi - 1)/Qphi;
[Theta, Phi] = ndgrid(theta, phi);
rho = problem.b*sin(Theta);
X = rho.*cos(Phi);
Y = rho.*sin(Phi);
Z = problem.a*cos(Theta);
g = sph_coords(problem.type, X(:)/problem.semifocal, ...
    Y(:)/problem.semifocal, Z(:)/problem.semifocal);
[E, H] = incident_plane_wave(problem.km, problem.n_m, ...
    problem.theta_inc, problem.phi_inc, problem.pol, X, Y, Z);
Eeta = reshape(sum(E.*g.e_eta, 2), size(Theta));
Ephi = reshape(sum(E.*g.e_phi, 2), size(Theta));
Heta = reshape(sum(H.*g.e_eta, 2), size(Theta));
Hphi = reshape(sum(H.*g.e_phi, 2), size(Theta));
incident = struct('Eeta', fft(Eeta, [], 2)/Qphi, ...
    'Ephi', fft(Ephi, [], 2)/Qphi, ...
    'Heta', fft(Heta, [], 2)/Qphi, ...
    'Hphi', fft(Hphi, [], 2)/Qphi, 'Qphi', Qphi);
end

function [g, theta, weights] = surface_meridian(problem, Qtheta)
[nodes, quadrature_weights] = gauss_legendre(Qtheta);
theta = (pi/2)*(nodes + 1);
weights = (pi/2)*quadrature_weights.*sin(theta);
eta = cos(theta);
if strcmp(problem.type, 'prolate')
    alpha = problem.xi0^2 - 1;
else
    alpha = problem.xi0^2 + 1;
end
rho = sqrt(alpha)*sqrt(max(0, 1 - eta.^2));
g = sph_coords(problem.type, rho, zeros(size(rho)), ...
    problem.xi0*eta);
end

function check = incident_tail_check(sol)
norms = [sol.modes.rhs_norm];
shell = abs([sol.modes.m]) >= max(0, sol.M - 2);
scale = max(norms);
if scale == 0
    value = 0;
else
    value = max(norms(shell))/scale;
end
check = make_check(true, value, value <= sol.opts.tol, ...
    struct('shell_norm', max(norms(shell)), 'maximum_norm', scale));
end

function check = qphi_check(sol, refined_Qphi)
[g, theta, weights] = surface_meridian(sol, sol.Qtheta);
incident = incident_spectrum(sol, theta, refined_Qphi);
maximum_rhs = max([sol.modes.rhs_norm]);
value = 0;
for mu = 0:sol.M
    data = sol.mu_data{mu + 1};
    if isempty(data)
        continue
    end
    Sext = sph_angular(data.eg_ext, g).W0;
    mode_indices = find(abs([sol.modes.m]) == mu);
    for k = mode_indices
        mode = sol.modes(k);
        bin = mod(mode.m, refined_Qphi) + 1;
        rhs = -[Sext*(weights.*incident.Eeta(:, bin)); ...
                Sext*(weights.*incident.Ephi(:, bin)); ...
                Sext*(weights.*incident.Heta(:, bin)); ...
                Sext*(weights.*incident.Hphi(:, bin))];
        value = max(value, ...
            norm(rhs - mode.rhs)/max(maximum_rhs, realmin));
    end
end
check = make_check(true, value, value <= sol.opts.tol, ...
    struct('refined_Qphi', refined_Qphi));
end

function check = modal_tail_check(sol)
count = sol.L + 10;
eta = cos((2*(1:count) - 1)*pi/(2*count)).';
if strcmp(sol.type, 'prolate')
    alpha = sol.xi0^2 - 1;
else
    alpha = sol.xi0^2 + 1;
end
rho = sqrt(alpha)*sqrt(max(0, 1 - eta.^2));
g = sph_coords(sol.type, rho, zeros(size(rho)), sol.xi0*eta);
value = 0;
for mu = 0:sol.M
    data = sol.mu_data{mu + 1};
    if isempty(data)
        continue
    end
    Wext = sph_angular(data.eg_ext, g);
    Wint = sph_angular(data.eg_int, g);
    Vext = sph_radial_ode_eval(data.rod_ext, g);
    [Vint, internal_ok] = sph_radial(data.eg_int, 1, g, sol.c_int);
    if ~internal_ok
        check = make_check(true, Inf, false, ...
            struct('reason', 'internal radial tail'));
        return
    end
    mode_indices = find(abs([sol.modes.m]) == mu);
    for k = mode_indices
        mode = sol.modes(k);
        [Mext, Next] = sph_vecwave(Wext, Vext, g, mode.m, sol.c_ext);
        [Mint, Nint] = sph_vecwave(Wint, Vint, g, mode.m, sol.c_int);
        top = max(1, numel(mode.a) - 2):numel(mode.a);
        ext_E = modal_norm(mode.a(top), mode.b(top), Next, Mext, top);
        ext_H = sol.n_m*modal_norm(mode.a(top), mode.b(top), ...
            Mext, Next, top);
        int_E = modal_norm(mode.c(top), mode.d(top), Nint, Mint, top);
        int_H = abs(sol.n_p)*modal_norm(mode.c(top), mode.d(top), ...
            Mint, Nint, top);
        value = max(value, max([ext_E, int_E]/sol.E0_norm));
        value = max(value, ...
            max([ext_H, int_H]/(sol.n_m*sol.E0_norm)));
    end
end
check = make_check(true, value, value <= sol.opts.tol, struct());
end

function value = modal_norm(first_coefficient, second_coefficient, ...
    first_wave, second_wave, rows)
eta = first_coefficient(:).'*first_wave.eta(rows, :) + ...
    second_coefficient(:).'*second_wave.eta(rows, :);
xi = first_coefficient(:).'*first_wave.xi(rows, :) + ...
    second_coefficient(:).'*second_wave.xi(rows, :);
phi = first_coefficient(:).'*first_wave.phi(rows, :) + ...
    second_coefficient(:).'*second_wave.phi(rows, :);
value = max(hypot(hypot(abs(eta), abs(xi)), abs(phi)));
end

function check = surface_residual_check(sol)
eta_count = sol.L + 10;
phi_count = 2*(2*sol.M + 5);
[coarse, coarse_detail] = boundary_residual(sol, eta_count, phi_count);
values = coarse;
grids = [eta_count, phi_count];
for refinement = 1:3
    refined_eta_count = 2*eta_count;
    refined_phi_count = 2*phi_count;
    [refined, refined_detail] = boundary_residual(sol, ...
        refined_eta_count, refined_phi_count);
    values(end + 1) = refined; %#ok<AGROW>
    grids(end + 1, :) = [refined_eta_count, refined_phi_count]; %#ok<AGROW>
    growth_change = max(0, refined - coarse);
    growth_limit = 0.1*max(coarse, 64*eps);
    growth_ok = refined <= 1.1*max(coarse, 64*eps);
    if growth_ok
        break
    end
    if refinement < 3
        eta_count = refined_eta_count;
        phi_count = refined_phi_count;
        coarse = refined;
        coarse_detail = refined_detail;
    end
end
passed = max(values) <= sol.opts.tol && growth_ok;
details = struct('coarse', coarse_detail, 'refined', refined_detail, ...
    'growth_ratio', refined/max(coarse, 64*eps), ...
    'growth_change', growth_change, ...
    'growth_change_limit', growth_limit, ...
    'coarse_grid', [eta_count, phi_count], ...
    'refined_grid', [refined_eta_count, refined_phi_count], ...
    'grid_history', grids, 'value_history', values);
check = make_check(true, max(values), passed, details);
end

function [value, detail] = boundary_residual(sol, eta_count, phi_count)
eta = cos((2*(1:eta_count) - 1)*pi/(2*eta_count)).';
phi = 2*pi*((0:phi_count - 1) + 0.5)/phi_count;
[Eta, Phi] = ndgrid(eta, phi);
sin_theta = sqrt(max(0, 1 - Eta.^2));
rho = sol.b*sin_theta;
X = rho.*cos(Phi);
Y = rho.*sin(Phi);
Z = sol.a*Eta;
[Eout, Hout] = spheroid_eval(sol, X, Y, Z, 'total', 'out');
[Ein, Hin] = spheroid_eval(sol, X, Y, Z, 'total', 'in');
normal_scale = min(sol.a, sol.b);
normal = [sin_theta(:).*cos(Phi(:))*(normal_scale/sol.b), ...
    sin_theta(:).*sin(Phi(:))*(normal_scale/sol.b), ...
    Eta(:)*(normal_scale/sol.a)];
normal = normal./row_norm(normal);
dE = Eout - Ein;
dH = Hout - Hin;
dE_tangent = dE - sum(dE.*normal, 2).*normal;
dH_tangent = dH - sum(dH.*normal, 2).*normal;
electric = max(row_norm(dE_tangent))/sol.E0_norm;
magnetic = max(row_norm(dH_tangent))/ ...
    (sol.n_m*sol.E0_norm);
normal_jump = sum((sol.n_m^2*Eout - sol.n_p^2*Ein).*normal, 2);
normal_value = max(abs(normal_jump))/ ...
    (sol.n_m^2*sol.E0_norm);
value = max([electric, magnetic, normal_value]);
detail = struct('tangential_E', electric, 'tangential_H', magnetic, ...
    'normal_n2E', normal_value, 'maximum', value);
end

function check = comparison_check(base, candidate, points, base_signature)
candidate_signature = solution_signature(candidate, points);
Escale = max(base.E0_norm, max(row_norm(base_signature.E)));
Hscale = max(base.n_m*base.E0_norm, max(row_norm(base_signature.H)));
electric = max(row_norm(candidate_signature.E - base_signature.E))/Escale;
magnetic = max(row_norm(candidate_signature.H - base_signature.H))/Hscale;
cross_section = abs(candidate_signature.Csca - base_signature.Csca)/ ...
    max(abs(base_signature.Csca), realmin);
value = max([electric, magnetic, cross_section]);
details = struct('electric', electric, 'magnetic', magnetic, ...
    'Csca', cross_section, 'candidate_L', candidate.L, ...
    'candidate_M', candidate.M, 'candidate_Qtheta', candidate.Qtheta, ...
    'candidate_Qphi', candidate.Qphi);
check = make_check(true, value, value <= base.opts.tol, details);
end

function signature = solution_signature(sol, points)
[E, H] = spheroid_eval(sol, points.X, points.Y, points.Z, ...
    'total', 'auto');
signature = struct('E', E, 'H', H, ...
    'Csca', scattering_cross_section(sol));
end

function value = scattering_cross_section(sol)
theta_count = 2*sol.L + 12;
phi_count = 2*(2*sol.M + 5);
[eta, weights] = gauss_legendre(theta_count);
phi = 2*pi*(0:phi_count - 1)/phi_count;
[Eta, Phi] = ndgrid(eta, phi);
length_scale = max(sol.a, sol.b);
radius_ratio = max([(1e7 + 1)*sol.semifocal/length_scale, ...
    1e4/(sol.km*length_scale), 2]);
radius = length_scale*radius_ratio;
sin_theta = sqrt(max(0, 1 - Eta.^2));
X = radius*sin_theta.*cos(Phi);
Y = radius*sin_theta.*sin(Phi);
Z = radius*Eta;
[E, ~] = spheroid_eval(sol, X, Y, Z, 'scattered', 'auto');
scaled_E = radius_ratio*(E/sol.E0_norm);
intensity = reshape(row_norm(scaled_E).^2, size(Eta));
value = (2*pi/phi_count)*sum(weights.*sum(intensity, 2));
end

function values = row_norm(array)
values = hypot(hypot(abs(array(:, 1)), abs(array(:, 2))), ...
    abs(array(:, 3)));
end

function points = deterministic_points(problem)
stream = RandStream('mt19937ar', 'Seed', 41723);
count = 25;
eta = -0.8 + 1.6*rand(stream, 2*count, 1);
phi = 2*pi*rand(stream, 2*count, 1);
if strcmp(problem.type, 'prolate')
    xi_in = 1 + (problem.xi0 - 1)*(0.1 + 0.8*rand(stream, count, 1));
    xi_out = problem.xi0*(1.05 + 0.5*rand(stream, count, 1));
    xi = [xi_in; xi_out];
    alpha = xi.^2 - 1;
else
    xi_in = problem.xi0*(0.1 + 0.8*rand(stream, count, 1));
    xi_out = problem.xi0 + (0.05 + 0.5*rand(stream, count, 1))* ...
        max(1, problem.xi0);
    xi = [xi_in; xi_out];
    alpha = xi.^2 + 1;
end
rho = sqrt(alpha.*(1 - eta.^2));
points = struct('X', problem.semifocal*rho.*cos(phi), ...
    'Y', problem.semifocal*rho.*sin(phi), ...
    'Z', problem.semifocal*xi.*eta);
end

function check = special_function_check(mu_data)
passed = true;
value = 0;
mu_data = mu_data(~cellfun(@isempty, mu_data));
details = repmat(struct('mu', 0, 'external_rmax_delta', NaN, ...
    'internal_rmax_delta', NaN, 'external_condition', NaN, ...
    'internal_condition', NaN, 'radial', struct()), numel(mu_data), 1);
for k = 1:numel(mu_data)
    data = mu_data{k};
    external_ok = data.eg_ext.cond < 1e8 && ...
        data.eg_ext.rmax_delta <= 3e-13 && data.eg_ext.rmax_tail <= 1e-14;
    internal_ok = data.eg_int.cond < 1e8 && ...
        data.eg_int.rmax_delta <= 3e-13 && data.eg_int.rmax_tail <= 1e-14;
    radial_ok = isfield(data.radial_check, 'ok') && ...
        islogical(data.radial_check.ok) && isscalar(data.radial_check.ok) && ...
        data.radial_check.ok;
    passed = passed && external_ok && internal_ok && radial_ok;
    value = max(value, max([data.eg_ext.rmax_delta, ...
        data.eg_int.rmax_delta, data.eg_ext.rmax_tail, ...
        data.eg_int.rmax_tail]));
    details(k) = struct('mu', data.mu, ...
        'external_rmax_delta', data.eg_ext.rmax_delta, ...
        'internal_rmax_delta', data.eg_int.rmax_delta, ...
        'external_condition', data.eg_ext.cond, ...
        'internal_condition', data.eg_int.cond, ...
        'radial', data.radial_check);
end
check = make_check(true, value, passed, details);
end

function checks = empty_checks()
blank = make_check(false, NaN, false, struct());
checks = struct('incident_tail', blank, 'qphi_independence', blank, ...
    'modal_tail', blank, 'surface_residual', blank, ...
    'qtheta_independence', blank, 'cutoff_stability', blank, ...
    'special_functions', blank);
end

function check = make_check(computed, value, passed, details)
check = struct('computed', logical(computed), 'value', value, ...
    'passed', logical(computed && passed), 'details', details);
end

function check = failed_check(exception)
details = struct('error_identifier', exception.identifier, ...
    'error_message', exception.message);
check = make_check(false, Inf, false, details);
end

function yes = all_checks_passed(checks)
names = fieldnames(checks);
yes = true;
for k = 1:numel(names)
    check = checks.(names{k});
    yes = yes && check.computed && check.passed;
end
end

function value = minimum_active_rank(modes)
ranks = [modes.numerical_rank];
ranks = ranks(isfinite(ranks));
if isempty(ranks)
    value = NaN;
else
    value = min(ranks);
end
end

function [M, L, Qtheta, Qphi, possible] = grow_M( ...
    M, L, Qtheta, Qphi, problem, opts)
candidate = M + ceil(0.25*M) + 2;
candidate_L = max(L, candidate);
possible = candidate <= opts.Mmax && candidate_L <= opts.Lmax;
if possible
    M = candidate;
    L = candidate_L;
    Qtheta = max(Qtheta, minimum_Qtheta(problem, L));
    Qphi = max(Qphi, minimum_Qphi(problem, M));
end
end

function [L, Qtheta, possible] = grow_L(L, Qtheta, problem, opts)
candidate = ceil(1.25*L) + 5;
possible = candidate <= opts.Lmax;
if possible
    L = candidate;
    Qtheta = max(Qtheta, minimum_Qtheta(problem, L));
end
end

function info = make_info(sol, checks, rounds, inside_envelope, ...
    completed, failure_reason, opts)
active_mask = ~cellfun(@isempty, sol.mu_data);
active_data = sol.mu_data(active_mask);
rmax_ext = nan(sol.M + 1, 1);
rmax_int = nan(sol.M + 1, 1);
xi_switch = nan(sol.M + 1, 1);
xi_start = nan(sol.M + 1, 1);
rmax_ext(active_mask) = cellfun(@(data) data.eg_ext.rmax, active_data);
rmax_int(active_mask) = cellfun(@(data) data.eg_int.rmax, active_data);
xi_switch(active_mask) = cellfun(@(data) data.rod_ext.xi_switch, active_data);
xi_start(active_mask) = cellfun(@(data) data.rod_ext.xi_start, active_data);
active = ~[sol.modes.analytically_zeroed];
if any(active)
    active_rcond = min([sol.modes(active).rcond]);
else
    active_rcond = Inf;
end
if checks.surface_residual.computed
    residual = checks.surface_residual.value;
else
    residual = NaN;
end
info = struct('route', 'spheroid', 'type', sol.type, 'L', sol.L, ...
    'M', sol.M, 'Qtheta', sol.Qtheta, 'Qphi', sol.Qphi, ...
    'rmax', max([rmax_ext(active_mask); rmax_int(active_mask)]), ...
    'rmax_by_mu', struct('external', rmax_ext, 'internal', rmax_int), ...
    'xi_switch', max(xi_switch(active_mask)), ...
    'xi_start', max(xi_start(active_mask)), ...
    'xi_switch_by_mu', xi_switch, 'xi_start_by_mu', xi_start, ...
    'residual', residual, 'rcond_min', active_rcond, ...
    'rcond_min_active', active_rcond, ...
    'solver_methods', {unique({sol.modes.solver_method})}, ...
    'numerical_rank_min', minimum_active_rank(sol.modes), ...
    'checks', checks, 'rounds', rounds, ...
    'special_function_diagnostics', checks.special_functions.details, ...
    'inside_envelope', inside_envelope, ...
    'validated', inside_envelope && completed && all_checks_passed(checks), ...
    'failure_reason', failure_reason, 'options', opts);
end

function yes = envelope_member(problem)
aspect = problem.a/problem.b;
n_relative = problem.n_p/problem.n_m;
loss_ratio = imag(n_relative)/real(n_relative);
scale = max(problem.a, problem.b);
x_ext = problem.km*scale;
yes = x_ext >= 0.02 && x_ext <= 250 && ...
    abs(problem.kp)*scale <= 250 && ...
    aspect >= 1/4 && aspect <= 4 && real(n_relative) >= 0.5 && ...
    real(n_relative) <= 2 && loss_ratio <= 0.1;
end
