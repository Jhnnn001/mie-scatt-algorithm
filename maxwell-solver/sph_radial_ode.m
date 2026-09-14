function [rod, chk] = sph_radial_ode(eg, c, xi0, xi_switch)
%SPH_RADIAL_ODE Build an arbitrary-scaled outgoing radial solution.

narginchk(4, 4);
validate_inputs(eg, c, xi0, xi_switch);
c = double(c);
xi0 = double(xi0);
xi_switch = max(double(xi_switch), xi0);
lambda = double(eg.lambda(:));
mu = double(eg.mu);
sigma = radial_sign(eg.type);
rel_tol = 1e-12;
abs_tol = 1e-12;

comparison_start = max([2 * xi_switch, 3, c / 3, ...
    sqrt(max(abs(lambda)))/c]);
coarse = integrate_from_start(eg.type, lambda, mu, c, sigma, ...
    xi0, comparison_start, rel_tol, abs_tol);
converged = false;
for attempt = 1:12
    xi_start = 1.5 * comparison_start;
    refined = integrate_from_start(eg.type, lambda, mu, c, sigma, ...
        xi0, xi_start, rel_tol, abs_tol);
    overlap_points = xi_switch + ...
        (comparison_start - xi_switch) * [1/4, 1/2, 3/4];
    check_points = [xi0, overlap_points];
    [log_error, field_error] = compare_solutions(coarse, refined, ...
        check_points, c);
    [ode_error, ode_floor] = solution_residual(refined, ...
        overlap_points, sigma);
    [asymptotic_error, asymptotic_floor] = asymptotic_residual( ...
        refined.coefficients, ...
        refined.stop, refined.argument_scale, lambda, mu, c, sigma, ...
        [xi_start, 1.5 * xi_start, 2 * xi_start]);
    asymptotic_floor = max(asymptotic_floor, ...
        32 * max(refined.tail));
    diagnostics_finite = all(isfinite([ode_error, ode_floor, ...
        asymptotic_error, asymptotic_floor, refined.tail(:).']));
    if log_error < 1e-10 && field_error < 1e-10 && diagnostics_finite
        converged = true;
        break
    end
    comparison_start = xi_start;
    coarse = refined;
end
if ~converged
    error('sph_radial_ode:StartConvergence', ...
        ['Outgoing start refinement failed (log %.3g, field %.3g, ' ...
        'ODE %.3g, asymptotic %.3g).'], log_error, field_error, ...
        ode_error, asymptotic_error);
end

backward_error_limit = 1e-9;
diagnostics_ok = diagnostics_finite && ...
    ode_error <= backward_error_limit && ...
    asymptotic_error <= backward_error_limit && ...
    ode_floor <= backward_error_limit && ...
    asymptotic_floor <= backward_error_limit;

rod = struct('type', eg.type, 'mu', mu, 'l', eg.l(:), ...
    'lambda', lambda, 'c', c, 'xi0', xi0, ...
    'xi_switch', xi_switch, 'comparison_start', comparison_start, ...
    'xi_start', xi_start, 'arbitrary_scale', true, ...
    'normalization', 'R(xi0)=1', 'normalization_xi', xi0, ...
    'solution', refined.solution, ...
    'qtilde_xi0', refined.qtilde_xi0, ...
    'asymptotic_coefficients', refined.coefficients, ...
    'asymptotic_stop', refined.stop, ...
    'asymptotic_argument_scale', refined.argument_scale, ...
    'asymptotic_U_start', refined.U_start, ...
    'RelTol', rel_tol, 'AbsTol', refined.abs_tol);
chk = struct('ok', diagnostics_ok, 'start_converged', true, ...
    'log_derivative_error', log_error, 'field_error', field_error, ...
    'backward_error_limit', backward_error_limit, ...
    'ode_residual', ode_error, 'ode_residual_floor', ode_floor, ...
    'asymptotic_residual', asymptotic_error, ...
    'asymptotic_residual_floor', asymptotic_floor, ...
    'asymptotic_tail', max(refined.tail), ...
    'overlap_points', overlap_points, 'xi_switch', xi_switch, ...
    'comparison_start', comparison_start, 'xi_start', xi_start, ...
    'start_attempts', attempt);
end

function validate_inputs(eg, c, xi0, xi_switch)
required = {'type', 'mu', 'lambda', 'l', 'c'};
if ~isstruct(eg) || ~all(isfield(eg, required))
    error('sph_radial_ode:InvalidEigenData', ...
        'eg must be a structure returned by sph_eigen.');
end
valid_modes = isnumeric(eg.lambda) && isvector(eg.lambda) && ...
    ~isempty(eg.lambda) && isreal(eg.lambda) && ...
    all(isfinite(eg.lambda(:))) && isnumeric(eg.l) && ...
    isvector(eg.l) && numel(eg.l) == numel(eg.lambda) && ...
    all(isfinite(eg.l(:)) & eg.l(:) == fix(eg.l(:)));
valid_mu = isnumeric(eg.mu) && isscalar(eg.mu) && isreal(eg.mu) && ...
    isfinite(eg.mu) && eg.mu >= 0 && eg.mu == fix(eg.mu);
if ~valid_modes || ~valid_mu || any(eg.l(:) < eg.mu)
    error('sph_radial_ode:InvalidEigenData', ...
        'eg contains invalid retained modes.');
end
if ~(isnumeric(c) && isscalar(c) && isreal(c) && isfinite(c) && c > 0)
    error('sph_radial_ode:InvalidExteriorParameter', ...
        'c must be a positive finite real scalar.');
end
if ~(isnumeric(xi0) && isscalar(xi0) && isreal(xi0) && isfinite(xi0))
    error('sph_radial_ode:InvalidXi0', 'xi0 must be a finite real scalar.');
end
if ~(isnumeric(xi_switch) && isscalar(xi_switch) && isreal(xi_switch) && ...
        isfinite(xi_switch) && xi_switch > 1)
    error('sph_radial_ode:InvalidSwitch', ...
        'xi_switch must be a finite real scalar greater than 1.');
end
if ~isnumeric(eg.c) || ~isscalar(eg.c) || ~isreal(eg.c) || ...
        ~isfinite(eg.c) || abs(double(eg.c) - double(c)) > ...
        8 * eps * max(1, double(c))
    error('sph_radial_ode:ParameterMismatch', ...
        'eg must have been computed using the supplied exterior c.');
end
if strcmp(eg.type, 'prolate')
    if xi0 <= 1
        error('sph_radial_ode:InvalidXi0', ...
            'Prolate xi0 must be greater than 1.');
    end
elseif strcmp(eg.type, 'oblate')
    if xi0 <= 0
        error('sph_radial_ode:InvalidXi0', ...
            'Oblate xi0 must be positive.');
    end
else
    error('sph_radial_ode:InvalidEigenData', ...
        'eg.type must be prolate or oblate.');
end
end

function run = integrate_from_start(type, lambda, mu, c, sigma, ...
        xi0, xi_start, rel_tol, abs_tol)
[coefficients, stop, w0, tail, argument_scale] = outgoing_seed(lambda, mu, c, ...
    sigma, xi_start, 80);
states = [w0; complex(zeros(size(w0)))];
modes = numel(lambda);
state_abs_tol = [min(1, c)*abs_tol*ones(modes, 1); ...
    abs_tol*ones(modes, 1)];
options = odeset('RelTol', rel_tol, 'AbsTol', state_abs_tol);
solution = ode113(@(xi, state) log_rhs(xi, state, lambda, mu, ...
    c, sigma), [xi_start, xi0], states, options);
at_xi0 = deval(solution, xi0);
[U_start, ~, ~] = inverse_polynomial(coefficients, stop, xi_start, ...
    argument_scale);
run = struct('type', type, 'lambda', lambda, 'mu', mu, 'c', c, ...
    'sigma', sigma, 'xi0', xi0, 'xi_start', xi_start, ...
    'solution', solution, ...
    'qtilde_xi0', at_xi0(modes + 1:end), ...
    'coefficients', coefficients, 'stop', stop, ...
    'argument_scale', argument_scale, 'U_start', U_start, ...
    'tail', tail, 'rel_tol', rel_tol, 'abs_tol', state_abs_tol);
end

function [coefficients, stop, w0, tail, argument_scale] = outgoing_seed(lambda, mu, ...
        c, sigma, xi, count)
modes = numel(lambda);
coefficients = complex(zeros(count + 1, modes));
stop = zeros(modes, 1);
w0 = complex(zeros(modes, 1));
tail = Inf(modes, 1);
s = 1i * c;
argument_scale = 1;
scaled_argument = c < 1;
if scaled_argument
    argument_scale = c;
end
for mode = 1:modes
    a = complex(zeros(count + 1, 1));
    a(1) = 1;
    for n = 0:count - 1
        previous1 = 0;
        previous2 = 0;
        if n >= 1
            previous1 = a(n);
        end
        if n >= 2
            previous2 = a(n - 1);
            powers = (1:floor(n / 2)).';
            indices = n - 2*powers;
            if scaled_argument
                parity_weight = sigma.^powers .* c.^(2*powers);
            else
                parity_weight = sigma.^powers;
            end
            parity_sum = sum(parity_weight .* a(indices + 1));
        else
            parity_sum = 0;
        end
        if scaled_argument
            numerator = (sigma * c^2 - lambda(mode) + n * (n + 1)) * ...
                a(n + 1) + 2 * sigma * 1i * n * c^2 * previous1 - ...
                sigma * n * (n - 1) * c^2 * previous2 - ...
                mu^2 * parity_sum;
            a(n + 2) = numerator / (2i * (n + 1));
        else
            numerator = (sigma * c^2 - lambda(mode) + n * (n + 1)) * ...
                a(n + 1) + 2 * sigma * s * n * previous1 - ...
                sigma * n * (n - 1) * previous2 - mu^2 * parity_sum;
            a(n + 2) = numerator / (2 * s * (n + 1));
        end
    end
    best_score = Inf;
    best_w = NaN;
    for used = 4:count + 1
        active = a(1:used);
        [U, Ux, Uxx, w] = inverse_polynomial(active, used, xi, ...
            argument_scale);
        y = s - 1 / xi + w;
        yp = Uxx / U - (Ux / U)^2;
        residual = stable_residual(y, yp, lambda(mode), mu, c, ...
            sigma, xi);
        powers = (argument_scale * xi).^-(0:used - 1).';
        terms = active .* powers;
        latest_tail = abs(terms(end)) / max(sum(abs(terms)), realmin);
        score = max(residual, latest_tail);
        if isfinite(score) && score < best_score
            best_score = score;
            stop(mode) = used;
            best_w = w;
            tail(mode) = latest_tail;
        end
    end
    if stop(mode) == 0 || ~isfinite(real(best_w)) || ...
            ~isfinite(imag(best_w))
        error('sph_radial_ode:AsymptoticSeed', ...
            'The outgoing asymptotic seed could not be formed.');
    end
    coefficients(:, mode) = a;
    w0(mode) = best_w;
end
end

function derivative = log_rhs(xi, state, lambda, mu, c, sigma)
modes = numel(lambda);
w = state(1:modes);
wp = stable_w_rhs(w, lambda, mu, c, sigma, xi);
derivative = [wp; w - 1 / xi];
end

function wp = stable_w_rhs(w, lambda, mu, c, sigma, xi)
alpha = xi^2 - sigma;
s = 1i * c;
wp = (lambda - sigma * c^2) / alpha - ...
    2 * sigma * s / (alpha * xi) + 2 * sigma / (alpha * xi^2) + ...
    sigma * mu^2 / alpha^2 - 2 * s * w - w.^2 - ...
    2 * sigma * w / (alpha * xi);
end

function [log_error, field_error] = compare_solutions(coarse, refined, ...
        xi, c)
[yc, Rc, Rpc] = normalized_values(coarse, xi);
[yr, Rr, Rpr] = normalized_values(refined, xi);
log_scale = max(abs(yr), c);
log_error = max(abs(yc - yr) ./ log_scale, [], 'all');
numerator = abs(Rc - Rr) + abs(Rpc - Rpr) / c;
denominator = max(abs(Rr) + abs(Rpr) / c, realmin);
field_error = max(numerator ./ denominator, [], 'all');
end

function [y, R, Rp] = normalized_values(run, xi)
state = deval(run.solution, xi);
modes = numel(run.lambda);
w = state(1:modes, :);
y = 1i * run.c - 1 ./ xi + w;
qtilde = state(modes + 1:end, :);
phase = exp(1i * run.c * (xi - run.xi0));
R = phase .* exp(qtilde - run.qtilde_xi0);
R(:, xi == run.xi0) = 1;
Rp = y .* R;
end

function [error_value, floor_value] = solution_residual(run, xi, sigma)
modes = numel(run.lambda);
error_value = 0;
floor_value = 32 * run.rel_tol;
for point = 1:numel(xi)
    available = min(xi(point) - run.xi0, run.xi_start - xi(point));
    derivative_step = max(5e-4, eps^(1/5) * max(1, abs(xi(point))));
    h = min([derivative_step, 0.01 / run.c, available / 4]);
    if ~(isfinite(h) && h > 0)
        error('sph_radial_ode:ResidualStep', ...
            'The radial residual requires a strictly interior point.');
    end
    offsets = [-2, -1, 1, 2];
    values = deval(run.solution, xi(point) + h * offsets);
    w = values(1:modes, :);
    derivative_h = (w(:, 1) - 8 * w(:, 2) + ...
        8 * w(:, 3) - w(:, 4)) / (12 * h);
    h = h / 2;
    values = deval(run.solution, xi(point) + h * offsets);
    w = values(1:modes, :);
    derivative_half = (w(:, 1) - 8 * w(:, 2) + ...
        8 * w(:, 3) - w(:, 4)) / (12 * h);
    wp = (16 * derivative_half - derivative_h) / 15;
    center = deval(run.solution, xi(point));
    y = 1i * run.c - 1 / xi(point) + center(1:modes);
    yp = 1 / xi(point)^2 + wp;
    [residual, scale] = radial_backward_residual(y, yp, run.lambda, ...
        run.mu, run.c, sigma, xi(point));
    if any(~isfinite(residual)) || any(~isfinite(scale))
        error('sph_radial_ode:NonfiniteResidual', ...
            'The radial ODE residual is nonfinite.');
    end
    error_value = max(error_value, max(residual));
    alpha = xi(point)^2 - sigma;
    estimate = alpha * abs(derivative_half - derivative_h) ./ ...
        (15 * scale);
    if any(~isfinite(estimate))
        error('sph_radial_ode:NonfiniteResidual', ...
            'The radial ODE residual floor is nonfinite.');
    end
    floor_value = max(floor_value, max(estimate));
end
end

function [error_value, floor_value] = asymptotic_residual( ...
        coefficients, stop, argument_scale, lambda, ...
        mu, c, sigma, xi)
[y, yp] = asymptotic_data(coefficients, stop, argument_scale, c, xi);
error_value = 0;
for point = 1:numel(xi)
    value = stable_residual(y(:, point), yp(:, point), lambda, mu, ...
        c, sigma, xi(point));
    if any(~isfinite(value))
        error('sph_radial_ode:NonfiniteResidual', ...
            'The radial asymptotic residual is nonfinite.');
    end
    error_value = max(error_value, max(value));
end
floor_value = 32 * eps;
end

function value = stable_residual(y, yp, lambda, mu, c, sigma, xi)
[value, ~] = radial_backward_residual(y, yp, lambda, mu, c, sigma, xi);
end

function [value, scale] = radial_backward_residual( ...
        y, yp, lambda, mu, c, sigma, xi)
alpha = xi^2 - sigma;
terms = [abs(alpha * yp), abs(alpha * y.^2), abs(2 * xi * y), ...
    repmat(abs(c^2 * xi^2), size(y)), abs(lambda), ...
    repmat(abs(mu^2 / alpha), size(y))];
scale = max(sum(terms, 2), realmin);
numerator = abs(alpha * (yp + y.^2) + 2 * xi * y + ...
    c^2 * xi^2 - lambda - sigma * mu^2 / alpha);
value = numerator ./ scale;
end

function [y, yp] = asymptotic_data(coefficients, stop, ...
        argument_scale, c, xi)
[U, Ux, Uxx] = inverse_polynomial(coefficients, stop, xi, ...
    argument_scale);
y = 1i * c + Ux ./ U;
yp = Uxx ./ U - (Ux ./ U).^2;
end

function [U, Ux, Uxx, w] = inverse_polynomial( ...
        coefficients, stop, xi, argument_scale)
xi = xi(:).';
modes = size(coefficients, 2);
U = complex(zeros(modes, numel(xi)));
Ux = U;
Uxx = U;
w = U;
for mode = 1:modes
    used = stop(min(mode, numel(stop)));
    if isscalar(stop)
        used = stop;
    end
    a = coefficients(1:used, min(mode, size(coefficients, 2)));
    n = (0:used - 1).';
    q = 1 ./ (argument_scale * xi);
    P = a.' * q.^n;
    if used >= 2
        Pq = (n(2:end) .* a(2:end)).' * q.^(n(2:end) - 1);
    else
        Pq = zeros(size(P));
    end
    if used >= 3
        Pqq = (n(3:end) .* (n(3:end) - 1) .* a(3:end)).' * ...
            q.^(n(3:end) - 2);
    else
        Pqq = zeros(size(P));
    end
    U(mode, :) = q .* P;
    Ux(mode, :) = -argument_scale * ...
        (q.^2 .* P + q.^3 .* Pq);
    Uxx(mode, :) = argument_scale^2 * ...
        (2 * q.^3 .* P + 4 * q.^4 .* Pq + q.^5 .* Pqq);
    w(mode, :) = -argument_scale * q.^2 .* Pq ./ P;
end
end

function sigma = radial_sign(type)
if strcmp(type, 'prolate')
    sigma = 1;
else
    sigma = -1;
end
end
