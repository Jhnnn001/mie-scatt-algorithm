function [V, ok] = sph_radial(eg, kind, g, c)
%SPH_RADIAL Evaluate spheroidal radial functions and weighted derivatives.

narginchk(4, 4);
required = {'type', 'mu', 'c', 'l', 'lambda', 'd', 'rmax'};
if ~isstruct(eg) || ~isscalar(eg) || ~all(isfield(eg, required))
    error('sph_radial:InvalidEigenData', ...
        'eg must be a structure returned by sph_eigen.');
end
if ~(isnumeric(kind) && isscalar(kind) && ismember(kind, 1:3))
    error('sph_radial:InvalidKind', 'kind must be 1, 2, or 3.');
end
if ~isstruct(g) || ~all(isfield(g, {'xi', 'alpha'}))
    error('sph_radial:InvalidCoordinates', 'g must contain xi and alpha.');
end
validateattributes(c, {'numeric'}, {'finite', 'scalar'}, mfilename, 'c');
c = double(c);
if ~(isnumeric(eg.c) && isscalar(eg.c) && ...
        isfinite(real(eg.c)) && isfinite(imag(eg.c))) || ...
        abs(eg.c - c) > 16*eps*max([1, abs(eg.c), abs(c)])
    error('sph_radial:ParameterMismatch', ...
        'c must match the parameter used to compute eg.');
end
validate_eigen_data(eg);
if kind ~= 1 && (~isreal(c) || c <= 0)
    error('sph_radial:InvalidExteriorParameter', ...
        'Kinds 2 and 3 require a positive real c.');
end
xi = g.xi(:).';
alpha = g.alpha(:).';
if numel(xi) ~= numel(alpha) || ~isreal(xi) || ~isreal(alpha) || ...
        any(~isfinite(xi)) || any(~isfinite(alpha))
    error('sph_radial:InvalidCoordinates', ...
        'xi and alpha must be equally sized finite real values.');
end
if strcmp(eg.type, 'prolate')
    if any(xi < 1) || any(alpha < 0)
        error('sph_radial:InvalidCoordinates', ...
            'Prolate coordinates require xi >= 1 and alpha >= 0.');
    end
elseif strcmp(eg.type, 'oblate')
    if any(xi < 0) || any(alpha < 1)
        error('sph_radial:InvalidCoordinates', ...
            'Oblate coordinates require xi >= 0 and alpha >= 1.');
    end
else
    error('sph_radial:InvalidEigenData', ...
        'eg.type must be prolate or oblate.');
end
if kind ~= 1 && any(xi <= 1)
    error('sph_radial:InvalidSecondKindDomain', ...
        'Kinds 2 and 3 require xi > 1.');
end

[coefficient, ok, well_conditioned] = series_coefficients(eg);
modes = size(eg.d, 2);
V = empty_output(modes, numel(xi));
if ~all(well_conditioned)
    if kind == 1
        V = regular_ivp(eg, g, c);
        ok = all(isfinite_output(V));
    else
        V = invalid_output(modes, numel(xi));
        ok = false;
    end
    return
end
if kind == 1 && strcmp(eg.type, 'oblate') && any(xi < 1)
    inner = xi < 1;
    [H, Hp, Hpp] = inner_sums(coefficient, eg.mu, c, xi(inner));
    Vinner = oblate_inner(H, Hp, Hpp, eg.mu, xi(inner), alpha(inner));
    V = assign_points(V, Vinner, inner);
end

outer = strcmp(eg.type, 'prolate') | xi >= 1;
if any(outer)
    bessel_kind = {'j', 'y', 'h1'};
    orders = eg.mu + (0:eg.rmax).';
    [z, dz, d2z] = sph_bessel(bessel_kind{kind}, orders, c*xi(outer));
    basis1 = c*dz;
    basis2 = c^2*d2z;
    G = coefficient.'*z;
    Gp = coefficient.'*basis1;
    Gpp = coefficient.'*basis2;
    Vouter = outer_weighted(G, Gp, Gpp, eg.mu, xi(outer), ...
        alpha(outer), eg.type);
    V = assign_points(V, Vouter, outer);
    if kind ~= 1
        ok = ok && tail_converged(coefficient, z, basis1, basis2, ...
            G, Gp, Gpp, eg.mu, eg.l);
    end
end
ok = ok && all(isfinite_output(V));
end

function validate_eigen_data(eg)
valid_indices = isnumeric(eg.mu) && isreal(eg.mu) && isscalar(eg.mu) && ...
    isfinite(eg.mu) && eg.mu >= 0 && eg.mu == fix(eg.mu) && ...
    isnumeric(eg.rmax) && isreal(eg.rmax) && isscalar(eg.rmax) && ...
    isfinite(eg.rmax) && eg.rmax >= 0 && eg.rmax == fix(eg.rmax);
if ~valid_indices || ~(isnumeric(eg.d) && ismatrix(eg.d) && ...
        ~isempty(eg.d) && size(eg.d, 1) == eg.rmax + 1)
    error('sph_radial:InvalidEigenData', ...
        'eg contains invalid modal indices or coefficient dimensions.');
end
modes = size(eg.d, 2);
valid_l = isnumeric(eg.l) && isreal(eg.l) && iscolumn(eg.l) && ...
    numel(eg.l) == modes && all(isfinite(eg.l)) && ...
    all(eg.l == fix(eg.l)) && all(eg.l >= eg.mu);
valid_lambda = isnumeric(eg.lambda) && iscolumn(eg.lambda) && ...
    numel(eg.lambda) == modes && ...
    all(isfinite(real(eg.lambda))) && all(isfinite(imag(eg.lambda)));
finite_d = all(isfinite(real(eg.d(:)))) && all(isfinite(imag(eg.d(:))));
if ~(valid_l && valid_lambda && finite_d)
    error('sph_radial:InvalidEigenData', ...
        'eg modal data must be finite and dimensionally consistent.');
end
end

function [coefficient, ok, well_conditioned] = series_coefficients(eg)
r = (0:eg.rmax).';
if size(eg.d, 1) ~= numel(r) || size(eg.d, 2) ~= numel(eg.l)
    error('sph_radial:InvalidEigenData', ...
        'eg.d dimensions must match eg.rmax and eg.l.');
end
log_weight = 0.5*(log((2*eg.mu + 2*r + 1)/(2*eg.mu + 1)) + ...
    gammaln(r + 2*eg.mu + 1) - gammaln(r + 1) - ...
    gammaln(2*eg.mu + 1));
weight = exp(log_weight - max(log_weight));
weighted = eg.d.*weight;
denominator = sum(weighted, 1);
denominator_scale = sum(abs(weighted), 1);
stable = abs(denominator) > eps*denominator_scale;
well_conditioned = stable & denominator_scale <= 1e4*abs(denominator);
safe_denominator = denominator;
safe_denominator(~stable) = 1;
exponent = r + eg.mu - eg.l(:).';
roots_i = [1, 1i, -1, -1i];
phase = reshape(roots_i(mod(exponent(:), 4) + 1), size(exponent));
coefficient = weighted.*phase./safe_denominator;
ok = all(stable) && all(isfinite(real(coefficient(:)))) && ...
    all(isfinite(imag(coefficient(:))));
end

function V = invalid_output(modes, points)
values = complex(nan(modes, points));
V = struct('V0', values, 'V1', values, 'V2', values, ...
    'Va', values, 'Vm', values);
end

function V = regular_ivp(eg, g, c)
xi = g.xi(:).';
alpha = g.alpha(:).';
if strcmp(eg.type, 'prolate')
    [R, Rp] = prolate_regular(eg, c, xi, alpha);
    coefficient = eg.lambda - c^2*xi.^2 + eg.mu^2./alpha;
else
    [R, Rp] = oblate_regular(eg, c, xi);
    coefficient = eg.lambda - c^2*xi.^2 - eg.mu^2./alpha;
end
Rpp_weighted = sqrt(alpha).*(-2*xi.*Rp + coefficient.*R);
V.V0 = R;
V.V1 = sqrt(alpha).*Rp;
V.V2 = Rpp_weighted;
if eg.mu == 0
    V.Va = complex(zeros(size(R)));
    V.Vm = V.Va;
else
    V.Va = R./sqrt(alpha);
    V.Vm = eg.mu*V.Va;
end
if strcmp(eg.type, 'prolate') && any(xi == 1)
    endpoint = xi == 1;
    V.V0(:, endpoint) = double(eg.mu == 0);
    V.V1(:, endpoint) = double(eg.mu == 1);
    V.V2(:, endpoint) = -double(eg.mu == 1);
    V.Va(:, endpoint) = double(eg.mu == 1);
    V.Vm(:, endpoint) = eg.mu*V.Va(:, endpoint);
end
end

function [R, Rp] = prolate_regular(eg, c, xi, alpha)
mu = eg.mu;
delta = min(1e-4, ...
    0.05/(1 + abs(c)^2 + max(abs(eg.lambda))));
if delta <= 64*eps
    error('sph_radial:RegularIVPUnavailable', ...
        'The prolate regular IVP start is below floating-point resolution.');
end
[F0, Fp0] = prolate_series(eg, c, delta);
integration_end = max(2, max(xi));
options = odeset('RelTol', 3e-13, 'AbsTol', 3e-13);
solution = ode113(@(x, y) prolate_rhs(x, y, eg.lambda, mu, c), ...
    [1 + delta, integration_end], [F0; Fp0], options);
F = complex(zeros(numel(eg.lambda), numel(xi)));
Fp = F;
near = xi <= 1 + delta;
if any(near)
    [F(:, near), Fp(:, near)] = prolate_series(eg, c, xi(near) - 1);
end
if any(~near)
    y = deval(solution, xi(~near));
    modes = numel(eg.lambda);
    F(:, ~near) = y(1:modes, :);
    Fp(:, ~near) = y(modes + 1:end, :);
end
power = alpha.^(mu/2);
R = power.*F;
Rp = complex(zeros(size(R)));
positive = alpha > 0;
Rp(:, positive) = power(positive).*(Fp(:, positive) + ...
    mu*xi(positive)./alpha(positive).*F(:, positive));
end

function [F, Fp] = prolate_series(eg, c, t)
count = 30;
modes = numel(eg.lambda);
a = complex(zeros(count + 1, modes));
a(1, :) = 1;
q0 = c^2 + eg.mu*(eg.mu + 1) - eg.lambda(:).';
for k = 0:count - 1
    previous = complex(zeros(1, modes));
    second_previous = previous;
    if k >= 1
        previous = a(k, :);
    end
    if k >= 2
        second_previous = a(k - 1, :);
    end
    numerator = (k*(k + 2*eg.mu + 1) + q0).*a(k + 1, :) + ...
        2*c^2*previous + c^2*second_previous;
    a(k + 2, :) = -numerator/(2*(k + 1)*(k + eg.mu + 1));
end
t = t(:).';
powers = t.^((0:count).');
F = a.'*powers;
orders = (1:count).';
Fp = (orders.*a(2:end, :)).'*t.^((0:count - 1).');
end

function dy = prolate_rhs(xi, y, lambda, mu, c)
modes = numel(lambda);
F = y(1:modes);
Fp = y(modes + 1:end);
alpha = xi^2 - 1;
Fpp = (-2*(mu + 1)*xi*Fp - ...
    (c^2*xi^2 + mu*(mu + 1) - lambda).*F)/alpha;
dy = [Fp; Fpp];
end

function [R, Rp] = oblate_regular(eg, c, xi)
modes = numel(eg.lambda);
even = mod(eg.l - eg.mu, 2) == 0;
R0 = double(even);
Rp0 = double(~even);
integration_end = max(1, max(xi));
options = odeset('RelTol', 3e-13, 'AbsTol', 3e-13);
solution = ode113(@(x, y) oblate_rhs(x, y, eg.lambda, eg.mu, c), ...
    [0, integration_end], [R0; Rp0], options);
y = deval(solution, xi);
R = y(1:modes, :);
Rp = y(modes + 1:end, :);
end

function dy = oblate_rhs(xi, y, lambda, mu, c)
modes = numel(lambda);
R = y(1:modes);
Rp = y(modes + 1:end);
alpha = xi^2 + 1;
Rpp = (-2*xi*Rp + ...
    (lambda - c^2*xi^2 - mu^2/alpha).*R)/alpha;
dy = [Rp; Rpp];
end

function V = empty_output(modes, points)
values = complex(zeros(modes, points));
V = struct('V0', values, 'V1', values, 'V2', values, ...
    'Va', values, 'Vm', values);
end

function V = assign_points(V, values, points)
names = fieldnames(V);
for k = 1:numel(names)
    V.(names{k})(:, points) = values.(names{k});
end
end

function V = outer_weighted(G, Gp, Gpp, mu, xi, alpha, type)
if strcmp(type, 'prolate')
    if mu == 0
        V.V0 = G;
        V.V1 = sqrt(alpha).*Gp;
        V.V2 = alpha.^(3/2).*Gpp;
        V.Va = zeros(size(G));
        V.Vm = zeros(size(G));
        return
    end
    base = (sqrt(alpha)./xi).^(mu - 1)./xi;
    V.V0 = sqrt(alpha).*base.*G;
    V.V1 = alpha.*base.*Gp + mu*base.*G./xi;
    V.V2 = alpha.^2.*base.*Gpp + ...
        2*mu*alpha.*base.*Gp./xi + ...
        mu*(mu + 1 - 3*xi.^2).*base.*G./xi.^2;
    V.Va = base.*G;
    V.Vm = mu*V.Va;
    return
end

F = (alpha./xi.^2).^(mu/2);
L = -mu./(alpha.*xi);
Lprime = mu*(3*xi.^2 + 1)./(alpha.^2.*xi.^2);
V.V0 = F.*G;
V.V1 = sqrt(alpha).*F.*(Gp + L.*G);
V.V2 = alpha.^(3/2).*F.*(Gpp + 2*L.*Gp + ...
    (L.^2 + Lprime).*G);
if mu == 0
    V.Va = zeros(size(G));
    V.Vm = zeros(size(G));
else
    V.Va = V.V0./sqrt(alpha);
    V.Vm = mu*V.Va;
end
end

function V = oblate_inner(H, Hp, Hpp, mu, xi, alpha)
F = alpha.^(mu/2);
L = mu*xi./alpha;
Lprime = mu*(1 - xi.^2)./alpha.^2;
V.V0 = F.*H;
V.V1 = sqrt(alpha).*F.*(Hp + L.*H);
V.V2 = alpha.^(3/2).*F.*(Hpp + 2*L.*Hp + ...
    (L.^2 + Lprime).*H);
if mu == 0
    V.Va = zeros(size(H));
    V.Vm = zeros(size(H));
else
    V.Va = V.V0./sqrt(alpha);
    V.Vm = mu*V.Va;
end
end

function [H, Hp, Hpp] = inner_sums(coefficient, mu, c, xi)
r = (0:size(coefficient, 1) - 1).';
[basis0, basis1, basis2] = combined_basis(mu, r, c, xi);
H = coefficient.'*basis0;
Hp = coefficient.'*basis1;
Hpp = coefficient.'*basis2;
end

function [H, Hp, Hpp] = combined_basis(mu, r, c, xi)
orders = mu + r;
if mu == 0
    [H, dj, d2j] = sph_bessel('j', orders, c*xi);
    Hp = c*dj;
    Hpp = c^2*d2j;
    return
end

[T, dT, d2T] = sph_bessel_T(mu, r, c*xi);
if c == 0
    H = complex(zeros(size(T)));
    Hp = H;
    Hpp = H;
else
    log_c = log(c);
    H = exp(log(T) + mu*log_c);
    Hp = exp(log(dT) + (mu + 1)*log_c);
    Hpp = exp(log(d2T) + (mu + 2)*log_c);
end
bad = T == 0 | dT == 0 | d2T == 0 | ...
    invalid(H) | invalid(Hp) | invalid(Hpp);
[rows, columns] = find(bad);
for k = 1:numel(rows)
    [H(rows(k), columns(k)), Hp(rows(k), columns(k)), ...
        Hpp(rows(k), columns(k))] = combined_series(mu, r(rows(k)), ...
        c, xi(columns(k)));
end
end

function bad = invalid(value)
bad = ~isfinite(real(value)) | ~isfinite(imag(value));
end

function [value, derivative, second] = combined_series(mu, r, c, xi)
n = mu + r;
log_coefficient = 0.5*log(pi) - (n + 1)*log(2) - gammaln(n + 1.5);
if xi == 0
    value = double(r == 0)*scaled_power(log_coefficient, c, n);
    derivative = double(r == 1)*scaled_power(log_coefficient, c, n);
    second = double(r == 2)*2*scaled_power(log_coefficient, c, n);
    if r == 0
        second = -scaled_power(log_coefficient, c, n + 2)/(2*n + 3);
    end
    return
end

if c == 0
    term = double(n == 0)*exp(log_coefficient)*xi^r;
else
    term = exp(log_coefficient + n*log(c) + r*log(xi));
end
value = term;
derivative = r*term/xi;
second = r*(r - 1)*term/xi^2;
scales = [abs(value), abs(derivative), abs(second)];
consecutive = 0;
x = c*xi;
for k = 1:max(100, ceil(2*abs(x) + 50))
    term = term*(-x^2)/(2*k*(2*n + 2*k + 1));
    power = r + 2*k;
    derivative_term = power*term/xi;
    second_term = power*(power - 1)*term/xi^2;
    value = value + term;
    derivative = derivative + derivative_term;
    second = second + second_term;
    scales = max(scales, [abs(term), abs(derivative_term), abs(second_term)]);
    latest = [abs(term), abs(derivative_term), abs(second_term)];
    if all(latest <= 8*eps*max(scales, realmin))
        consecutive = consecutive + 1;
        if consecutive == 3
            return
        end
    else
        consecutive = 0;
    end
end
error('sph_radial:SeriesNoConvergence', ...
    'Combined scaled spherical-Bessel series did not converge.');
end

function value = scaled_power(log_coefficient, base, exponent)
if base == 0
    value = double(exponent == 0)*exp(log_coefficient);
else
    value = exp(log_coefficient + exponent*log(base));
end
end

function ok = tail_converged(coefficient, B0, B1, B2, S0, S1, S2, mu, l)
r = (0:size(coefficient, 1) - 1).';
ok = true;
for mode = 1:size(coefficient, 2)
    active = find(mod(r, 2) == mod(l(mode) - mu, 2));
    last = active(max(1, numel(active) - 4):end);
    scale = abs(coefficient(last, mode));
    tail0 = max(scale.*abs(B0(last, :)), [], 1);
    tail1 = max(scale.*abs(B1(last, :)), [], 1);
    tail2 = max(scale.*abs(B2(last, :)), [], 1);
    ok = ok && all(tail0 <= 1e-14*max(abs(S0(mode, :)), realmin)) && ...
        all(tail1 <= 1e-14*max(abs(S1(mode, :)), realmin)) && ...
        all(tail2 <= 1e-14*max(abs(S2(mode, :)), realmin));
end
ok = ok && all(isfinite(real(B0(:)))) && all(isfinite(imag(B0(:)))) && ...
    all(isfinite(real(B1(:)))) && all(isfinite(imag(B1(:)))) && ...
    all(isfinite(real(B2(:)))) && all(isfinite(imag(B2(:))));
end

function values = isfinite_output(V)
values = [isfinite(V.V0(:)); isfinite(V.V1(:)); ...
    isfinite(V.V2(:)); isfinite(V.Va(:)); isfinite(V.Vm(:))];
end
