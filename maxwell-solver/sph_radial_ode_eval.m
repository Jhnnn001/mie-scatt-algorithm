function V = sph_radial_ode_eval(rod, g)
%SPH_RADIAL_ODE_EVAL Evaluate an arbitrary-scaled outgoing radial basis.

narginchk(2, 2);
required = {'type', 'mu', 'lambda', 'c', 'xi0', 'xi_start', ...
    'arbitrary_scale', 'normalization_xi', 'solution', 'qtilde_xi0', ...
    'asymptotic_coefficients', 'asymptotic_stop', ...
    'asymptotic_argument_scale', ...
    'asymptotic_U_start'};
if ~isstruct(rod) || ~all(isfield(rod, required)) || ...
        ~isequal(rod.arbitrary_scale, true) || ...
        rod.normalization_xi ~= rod.xi0
    error('sph_radial_ode_eval:InvalidSolution', ...
        'rod must be an arbitrary-scaled structure from sph_radial_ode.');
end
if ~isstruct(g) || ~all(isfield(g, {'xi', 'alpha'}))
    error('sph_radial_ode_eval:InvalidCoordinates', ...
        'g must contain xi and alpha.');
end
xi = g.xi(:).';
alpha = g.alpha(:).';
if numel(xi) ~= numel(alpha) || ~isreal(xi) || ~isreal(alpha) || ...
        any(~isfinite(xi)) || any(~isfinite(alpha))
    error('sph_radial_ode_eval:InvalidCoordinates', ...
        'xi and alpha must be equally sized finite real values.');
end
if strcmp(rod.type, 'prolate')
    valid_geometry = xi > 1 & alpha > 0;
elseif strcmp(rod.type, 'oblate')
    valid_geometry = xi >= 0 & alpha >= 1;
else
    error('sph_radial_ode_eval:InvalidSolution', ...
        'rod.type must be prolate or oblate.');
end
if ~all(valid_geometry)
    error('sph_radial_ode_eval:InvalidCoordinates', ...
        'Coordinates are outside the radial domain.');
end
tolerance = 32 * eps * max(1, rod.xi_start);
if any(xi < rod.xi0 - tolerance)
    error('sph_radial_ode_eval:OutsideInterval', ...
        'xi must be greater than or equal to rod.xi0.');
end
xi = max(xi, rod.xi0);

modes = numel(rod.lambda);
valid_solution = numel(rod.qtilde_xi0) == modes && ...
    isnumeric(rod.asymptotic_argument_scale) && ...
    isscalar(rod.asymptotic_argument_scale) && ...
    isreal(rod.asymptotic_argument_scale) && ...
    isfinite(rod.asymptotic_argument_scale) && ...
    rod.asymptotic_argument_scale > 0 && ...
    isequal(size(rod.asymptotic_coefficients, 2), modes) && ...
    numel(rod.asymptotic_stop) == modes && ...
    numel(rod.asymptotic_U_start) == modes && ...
    all(rod.asymptotic_stop(:) >= 1) && ...
    all(rod.asymptotic_stop(:) == fix(rod.asymptotic_stop(:))) && ...
    all(rod.asymptotic_stop(:) <= ...
    size(rod.asymptotic_coefficients, 1));
if ~valid_solution
    error('sph_radial_ode_eval:InvalidSolution', ...
        'rod contains inconsistent mode or asymptotic data.');
end
stored_values = [rod.qtilde_xi0(:); ...
    rod.asymptotic_coefficients(:); rod.asymptotic_U_start(:)];
if any(~isfinite(real(stored_values))) || ...
        any(~isfinite(imag(stored_values)))
    error('sph_radial_ode_eval:NonfiniteOutput', ...
        'The stored radial ODE data contain a nonfinite value.');
end
R = complex(zeros(modes, numel(xi)));
Rp = R;
inside = xi <= rod.xi_start + tolerance;
if any(inside)
    bounded_xi = min(xi(inside), rod.xi_start);
    [sorted_xi, order] = sort(bounded_xi);
    state_sorted = deval(rod.solution, sorted_xi);
    state = zeros(size(state_sorted), 'like', state_sorted);
    state(:, order) = state_sorted;
    w = state(1:modes, :);
    y = 1i * rod.c - 1 ./ bounded_xi + w;
    qtilde = state(modes + 1:end, :);
    phase = exp(1i * rod.c * (bounded_xi - rod.xi0));
    R(:, inside) = phase .* exp(qtilde - rod.qtilde_xi0(:));
    R(:, inside & xi == rod.xi0) = 1;
    Rp(:, inside) = y .* R(:, inside);
end
outside = ~inside;
if any(outside)
    [U, Ux] = inverse_polynomial(rod.asymptotic_coefficients, ...
        rod.asymptotic_stop, xi(outside), ...
        rod.asymptotic_argument_scale);
    y = 1i * rod.c + Ux ./ U;
    start_phase = exp(1i * rod.c * (rod.xi_start - rod.xi0));
    R_start = start_phase .* exp(-rod.qtilde_xi0(:));
    phase = exp(1i * rod.c * (xi(outside) - rod.xi_start));
    ratio = U ./ rod.asymptotic_U_start(:);
    R(:, outside) = R_start .* phase .* ratio;
    Rp(:, outside) = y .* R(:, outside);
end

if strcmp(rod.type, 'prolate')
    Rpp = (-2 * xi .* Rp + ...
        (rod.lambda(:) - rod.c^2 * xi.^2 + ...
        rod.mu^2 ./ alpha) .* R) ./ alpha;
else
    Rpp = (-2 * xi .* Rp + ...
        (rod.lambda(:) - rod.c^2 * xi.^2 - ...
        rod.mu^2 ./ alpha) .* R) ./ alpha;
end
V0 = R;
V1 = sqrt(alpha) .* Rp;
V2 = alpha.^(3/2) .* Rpp;
if rod.mu == 0
    Va = complex(zeros(size(R)));
    Vm = Va;
else
    Va = R ./ sqrt(alpha);
    Vm = rod.mu * Va;
end
values = [V0(:); V1(:); V2(:); Va(:); Vm(:)];
if any(~isfinite(real(values))) || any(~isfinite(imag(values)))
    error('sph_radial_ode_eval:NonfiniteOutput', ...
        'The radial ODE evaluation produced a nonfinite value.');
end
V = struct('V0', V0, 'V1', V1, 'V2', V2, 'Va', Va, 'Vm', Vm);
end

function [U, Ux] = inverse_polynomial( ...
        coefficients, stop, xi, argument_scale)
xi = xi(:).';
modes = size(coefficients, 2);
U = complex(zeros(modes, numel(xi)));
Ux = U;
for mode = 1:modes
    used = stop(mode);
    a = coefficients(1:used, mode);
    n = (0:used - 1).';
    q = 1 ./ (argument_scale * xi);
    P = a.' * q.^n;
    if used >= 2
        Pq = (n(2:end) .* a(2:end)).' * q.^(n(2:end) - 1);
    else
        Pq = zeros(size(P));
    end
    U(mode, :) = q .* P;
    Ux(mode, :) = -argument_scale * ...
        (q.^2 .* P + q.^3 .* Pq);
end
end
