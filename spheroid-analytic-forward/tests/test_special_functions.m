function test_special_functions
% Spheroidal special-function checks.

test_zero_limit;
test_eigen_system;
test_roundoff_floor;
test_eigen_sensitivity;
test_angular_closed_forms;
test_angular_orthonormality;
test_angular_ode;
test_high_order_endpoints;
test_bessel_closed_forms_and_layout;
test_bessel_zero_values;
test_bessel_ode_and_wronskian;
test_scaled_bessel_identity;
test_scaled_bessel_zero_and_paths;
test_bessel_underflow_and_high_order;
test_radial_spherical_limit;
test_radial_ode;
test_outgoing_radial_ivp;
test_radial_weighted_endpoints;
test_radial_oblate_path_agreement;
test_radial_oblate_scaled_derivatives;
test_radial_high_mu_combined_scale;
test_radial_high_mu_far_point;
test_radial_high_c_regular_ivp;
test_radial_regular_ivp_scales;
test_radial_outgoing_tail_and_identity;
test_radial_tail_parity_and_normalization;
test_radial_outgoing_validation;
test_radial_validation_fails_closed;
test_validation;

fprintf('test_special_functions: PASS\n');
end

function test_eigen_sensitivity
must_error_id(@() sph_eigen('oblate', 0, 200 + 20i, 234), ...
    'sph_eigen:EigenSensitivity');
end

function test_zero_limit
mu = 2;
L = 7;
eg = sph_eigen('prolate', mu, 1e-6, L);
ell = (mu:L).';
assert(max(abs(eg.lambda - ell .* (ell + 1))) < 1e-8);
for k = 1:numel(ell)
    expected = zeros(size(eg.r));
    expected(ell(k) - mu + 1) = 1;
    assert(norm(eg.d(:, k) - expected) < 1e-10);
end
end

function test_eigen_system
type = 'oblate';
mu = 1;
c = 3 + 0.2i;
L = 6;
eg = sph_eigen(type, mu, c, L);

assert(isequal(eg.r, (0:eg.rmax).'));
assert(isequal(eg.l, (mu:L).'));
assert(isequal(size(eg.d), [eg.rmax + 1, L - mu + 1]));
assert(norm(eg.d.' * eg.d - eye(L - mu + 1), 'fro') < 2e-11);
assert(isfinite(eg.cond) && eg.cond < 1e8);
assert(abs(eg.cond / cond(eg.d) - 1) < 10 * eps);
assert(isfinite(eg.rmax_delta) && eg.rmax_delta <= 3e-13);
assert(isfinite(eg.rmax_tail) && eg.rmax_tail <= 1e-14);
assert(islogical(eg.rmax_floor_used) && isscalar(eg.rmax_floor_used));
initial_rmax = (L - mu) + ceil(abs(c)) + 40;
assert(eg.rmax >= ceil(1.5 * initial_rmax));

for k = 1:numel(eg.l)
    parity = mod(eg.l(k) - mu, 2);
    assert(all(eg.d(mod(eg.r, 2) ~= parity, k) == 0));
    target = eg.l(k) - mu + 1;
    assert(real(eg.d(target, k)) >= -1e-14);
end

A = angular_matrix(type, mu, c, eg.r);
residual = A * eg.d - eg.d .* eg.lambda.';
scale = max(1, norm(A, 'fro') * norm(eg.d, 'fro'));
assert(norm(residual, 'fro') / scale < 2e-13);
end

function test_roundoff_floor
c = 200 + 20i;
eg = sph_eigen('oblate', 0, c, 4);
initial_rmax = 4 + ceil(abs(c)) + 40;
assert(eg.rmax == ceil(1.5 * initial_rmax));
assert(eg.rmax_delta <= 3e-13);
assert(eg.rmax_tail <= 1e-14);
if eg.rmax_floor_used
    assert(eg.rmax_delta >= 1e-13);
else
    assert(eg.rmax_delta < 1e-13);
end
end

function test_angular_closed_forms
eta = [-1; -0.4; 0; 0.3; 1];
g = struct('eta', eta, 'beta', 1 - eta.^2);
W = sph_angular(sph_eigen('prolate', 0, 0, 2), g);

root_beta = sqrt(g.beta).';
eta_row = eta.';
expected_W0 = [sqrt(1/2) * ones(1, numel(eta)); ...
    sqrt(3/2) * eta_row; ...
    sqrt(5/8) * (3 * eta_row.^2 - 1)];
expected_W1 = [zeros(1, numel(eta)); ...
    sqrt(3/2) * root_beta; ...
    6 * sqrt(5/8) * eta_row .* root_beta];
expected_W2 = [zeros(2, numel(eta)); ...
    6 * sqrt(5/8) * root_beta.^3];
error_W0 = abs(W.W0 - expected_W0);
error_W1 = abs(W.W1 - expected_W1);
error_W2 = abs(W.W2 - expected_W2);
assert(max(error_W0(:)) < 2e-13);
assert(max(error_W1(:)) < 2e-13);
assert(max(error_W2(:)) < 2e-13);
assert(all(W.Wm(:) == 0));

g = struct('eta', [0.3; 1; -1], 'beta', [0.4; 0; 0]);
W = sph_angular(sph_eigen('prolate', 1, 0, 1), g);
seed = sqrt(3/4);
assert(max(abs(W.W0 - seed * sqrt(g.beta).')) < 2e-13);
assert(max(abs(W.W1 + seed * g.eta.')) < 2e-13);
assert(max(abs(W.W2 + seed * (g.eta.^2 + g.beta).')) < 2e-13);
assert(max(abs(W.Wm - seed)) < 2e-13);
assert(all(isfinite([W.W0(:); W.W1(:); W.W2(:); W.Wm(:)])));
end

function test_angular_orthonormality
for type = {'prolate', 'oblate'}
    for c = [0.5, 5, 50, 200]
        for mu = [0, 3, 40, 200]
            L = mu + 4;
            eg = sph_eigen(type{1}, mu, c, L);
            [eta, weight] = gauss_legendre(mu + eg.rmax + 3);
            W = sph_angular(eg, struct('eta', eta, ...
                'beta', 1 - eta.^2));
            gram = W.W0 * (weight .* W.W0.');
            assert(norm(gram - eye(L - mu + 1), 'fro') < 1e-10);
        end
    end
end

eg = sph_eigen('prolate', 2, 3 + 0.1i, 6);
[eta, weight] = gauss_legendre(eg.mu + eg.rmax + 3);
W = sph_angular(eg, struct('eta', eta, 'beta', 1 - eta.^2));
assert(norm(W.W0 * (weight .* W.W0.') - eye(5), 'fro') < 1e-10);
end

function test_angular_ode
for type = {'prolate', 'oblate'}
    mu = 2;
    c = 4 + 0.1i;
    eg = sph_eigen(type{1}, mu, c, 6);
    eta = linspace(-0.85, 0.85, 13).';
    beta = 1 - eta.^2;
    W = sph_angular(eg, struct('eta', eta, 'beta', beta));
    S = W.W0;
    Sp = W.W1 ./ sqrt(beta).';
    Spp = W.W2 ./ (beta.^(3/2)).';
    sigma = 1;
    if strcmp(type{1}, 'oblate')
        sigma = -1;
    end
    potential = eg.lambda + (-sigma * c^2 * eta.^2 - mu^2 ./ beta).';
    residual = beta.' .* Spp - 2 * eta.' .* Sp + potential .* S;
    scale = abs(beta.' .* Spp) + abs(2 * eta.' .* Sp) + ...
        abs(potential .* S);
    residual = abs(residual);
    assert(max(residual(:)) / max(1, max(scale(:))) < 1e-9);
end
end

function test_high_order_endpoints
mu = 200;
eg = sph_eigen('oblate', mu, 0.5, mu + 1);
eta = [-1; -0.2; 0; 0.8; 1];
W = sph_angular(eg, struct('eta', eta, 'beta', 1 - eta.^2));
assert(all(isfinite([W.W0(:); W.W1(:); W.W2(:); W.Wm(:)])));
end

function test_bessel_closed_forms_and_layout
x = [0.7, -0.4 + 0.3i];
[j, dj, d2j] = sph_bessel('j', [0; 1], x);
j0 = sin(x)./x;
j1 = sin(x)./x.^2 - cos(x)./x;
assert(isequal(size(j), [2, 2]));
assert_relative(j, [j0; j1], 2e-13);
assert_relative(dj(1, :), -j1, 2e-13);
assert_relative(d2j(1, :), ...
    -sin(x)./x - 2*cos(x)./x.^2 + 2*sin(x)./x.^3, 2e-12);

[y, ~, ~] = sph_bessel('y', [0; 1], x);
y0 = -cos(x)./x;
y1 = -cos(x)./x.^2 - sin(x)./x;
assert_relative(y, [y0; y1], 2e-13);
h = sph_bessel('h1', [0; 1], x);
assert_relative(h, j + 1i*y, 3e-13);
end

function test_bessel_zero_values
[j, dj, d2j] = sph_bessel('j', 0:4, 0);
assert(isequal(j, [1; 0; 0; 0; 0]));
assert(isequal(dj, [0; 1/3; 0; 0; 0]));
assert(isequal(d2j, [-1/3; 0; 2/15; 0; 0]));
must_error(@() sph_bessel('y', 0, 0));
must_error(@() sph_bessel('h1', 0, 0));
end

function test_bessel_ode_and_wronskian
n = [0; 1; 7; 30];
x = [0.8 + 0.2i, 3.2 - 0.4i];
[z, dz, d2z] = sph_bessel('j', n, x);
N = n(:);
X = x(:).';
terms = cat(3, X.^2.*d2z, 2*X.*dz, ...
    (X.^2 - N.*(N + 1)).*z);
residual = sum(terms, 3);
scale = sum(abs(terms), 3);
assert(max(abs(residual(:))./max(scale(:), realmin)) < 3e-12);

n = [0; 1; 8];
x = [0.7, 2.3, 250];
[j, dj] = sph_bessel('j', n, x);
[y, dy] = sph_bessel('y', n, x);
W = j.*dy - dj.*y;
assert_relative(W, repmat(1./x.^2, numel(n), 1), 2e-11);
end

function test_scaled_bessel_identity
mu = 4;
r = [0; 1; 4];
x = [0.2, 1.3 + 0.4i, -0.8 + 0.2i];
[T, dT, d2T] = sph_bessel_T(mu, r, x);
[j, dj, d2j] = sph_bessel('j', mu + r, x);
X = x(:).';
expected = j./X.^mu;
expected_d = dj./X.^mu - mu*j./X.^(mu + 1);
expected_d2 = d2j./X.^mu - 2*mu*dj./X.^(mu + 1) + ...
    mu*(mu + 1)*j./X.^(mu + 2);
assert(isequal(size(T), [numel(r), numel(x)]));
assert_relative(T, expected, 3e-12);
assert_relative(dT, expected_d, 3e-11);
assert_relative(d2T, expected_d2, 3e-10);
end

function test_scaled_bessel_zero_and_paths
mu = 3;
r = (0:3).';
[T, dT, d2T] = sph_bessel_T(mu, r, 0);
C = arrayfun(@(rr) inverse_odd_double_factorial(mu + rr), r);
assert_relative(T, [C(1); 0; 0; 0], 2e-15);
assert_relative(dT, [0; C(2); 0; 0], 2e-15);
assert_relative(d2T, [-C(1)/(2*mu + 3); 0; 2*C(3); 0], 2e-15);

x = 1e-8*[1, 1i, -1, -1i];
[T, dT] = sph_bessel_T(mu, [0; 1], x);
C0 = inverse_odd_double_factorial(mu);
C1 = inverse_odd_double_factorial(mu + 1);
assert_relative(T(1, :), C0*ones(size(x)), 2e-15);
assert_relative(T(2, :), C1*x, 2e-15);
assert_relative(dT(2, :), C1*ones(size(x)), 2e-15);

phase = exp(0.63i);
x = phase*[0.999e-3, 1.001e-3];
[T, dT] = sph_bessel_T(mu, [0; 2], x);
for k = 1:numel(x)
    [Ts, dTs] = scaled_bessel_series(mu, [0; 2], x(k));
    assert_relative(T(:, k), Ts, 2e-13);
    assert_relative(dT(:, k), dTs, 2e-13);
end
end

function test_bessel_underflow_and_high_order
mu = 100;
x = 0.002*exp(0.4i);
j = sph_bessel('j', mu, x);
assert(j == 0);
[T, dT, d2T] = sph_bessel_T(mu, 0, x);
[expected, expected_d, expected_d2] = scaled_bessel_series(mu, 0, x);
assert(T ~= 0);
assert_relative(T, expected, 3e-13);
assert_relative(dT, expected_d, 3e-13);
assert_relative(d2T, expected_d2, 3e-13);

n = 200;
x = 250;
j_nm1 = sph_bessel('j', n - 1, x);
j_n = sph_bessel('j', n, x);
j_np1 = sph_bessel('j', n + 1, x);
assert_relative(j_nm1 + j_np1, (2*n + 1)*j_n/x, 2e-12);
end

function test_radial_spherical_limit
mu = 0;
l = (0:4).';
c = 1e-6;
eg = unit_eigen('prolate', mu, l, 18, c);
xi = 3/c;
g = radial_geometry('prolate', xi);
[V1, ok1] = sph_radial(eg, 1, g, c);
[V2, ok2] = sph_radial(eg, 2, g, c);
j = sph_bessel('j', l, 3);
y = sph_bessel('y', l, 3);
assert(ok1 && ok2);
assert_relative(V1.V0, j, 2e-12);
assert_relative(V2.V0, y, 1e-8);
end

function test_radial_ode
mu = 1;
c = 1.2;
eg = sph_eigen('prolate', mu, c, 4);
xi = [1.1, 1.5, 2.3];
g = radial_geometry('prolate', xi);
[V, ok] = sph_radial(eg, 1, g, c);
assert(ok);
assert_radial_ode(eg, V, g, c, 'prolate', 2e-10);

mu = 2;
c = 0.7 + 0.2i;
eg = sph_eigen('oblate', mu, c, 4);
xi = [0, 0.2, 0.7, 1, 1.4];
g = radial_geometry('oblate', xi);
[V, ok] = sph_radial(eg, 1, g, c);
assert(ok);
assert(all(finite_radial_values(V)));
assert_radial_ode(eg, V, g, c, 'oblate', 2e-9);
end

function test_outgoing_radial_ivp
cases = { ...
    'prolate', 50, 0, 2, 1.03; ...
    'prolate', 200, 0, 2, 1.03; ...
    'prolate', 50, 40, 42, 1.03; ...
    'prolate', 200, 200, 200, 1.03; ...
    'oblate', 50, 0, 2, 0.577; ...
    'oblate', 200, 0, 2, 0.26; ...
    'oblate', 50, 40, 42, 0.26; ...
    'oblate', 200, 200, 200, 0.577};
rods = cell(size(cases, 1), 1);
eigen_data = rods;
for k = 1:size(cases, 1)
    type = cases{k, 1};
    c = cases{k, 2};
    mu = cases{k, 3};
    L = cases{k, 4};
    xi0 = cases{k, 5};
    eg = sph_eigen(type, mu, c, L);
    [rod, chk] = sph_radial_ode(eg, c, xi0, 1.5);
    assert(rod.arbitrary_scale && rod.normalization_xi == xi0);
    assert(strcmp(rod.normalization, 'R(xi0)=1'));
    assert(chk.ok && chk.start_converged && ...
        chk.log_derivative_error < 1e-10);
    assert(chk.field_error < 1e-10);
    backward_errors = [chk.ode_residual, chk.asymptotic_residual, ...
        chk.ode_residual_floor, chk.asymptotic_residual_floor];
    assert(chk.backward_error_limit <= 1e-9);
    assert(all(backward_errors <= 1e-9));
    V0 = sph_radial_ode_eval(rod, radial_geometry(type, xi0));
    assert(max(abs(V0.V0 - 1)) < 5e-11);
    assert(outgoing_fd_residual(rod, eg) < 2e-7);
    rods{k} = rod;
    eigen_data{k} = eg;
end

eg = sph_eigen('prolate', 0, 5, 2);
[~, chk] = sph_radial_ode(eg, 5, 2, 1.5);
overlap_lower_bound = max(2, chk.xi_switch);
assert(all(chk.overlap_points > overlap_lower_bound));
assert(all(chk.overlap_points < chk.comparison_start));

for k = 1:numel(rods)
    rod = rods{k};
    xi = 20 * rod.xi_start;
    V = sph_radial_ode_eval(rod, radial_geometry(rod.type, xi));
    assert(all(finite_radial_values(V)));
end

rod = rods{1};
eg = eigen_data{1};
xi = 20 * rod.xi_start;
g = radial_geometry(rod.type, xi);
V = sph_radial_ode_eval(rod, g);
y = (V.V1 / sqrt(g.alpha)) ./ V.V0;
sigma = 1;
A = (eg.lambda - sigma * rod.c^2) / (2i * rod.c);
B = -sigma / 2 - eg.lambda / (2 * rod.c^2);
y3 = 1i * rod.c - 1 / xi + A / xi^2 + B / xi^3;
assert(max(abs(y - y3)) * xi^3 < 0.2);

must_error_id(@() sph_radial_ode(eg, rod.c + 0.1i, ...
    rod.xi0, 1.5), 'sph_radial_ode:InvalidExteriorParameter');
must_error_id(@() sph_radial_ode(eigen_data{2}, 200, 0, 1.5), ...
    'sph_radial_ode:InvalidXi0');
must_error_id(@() sph_radial_ode_eval(rod, ...
    radial_geometry(rod.type, rod.xi0 - 1e-3)), ...
    'sph_radial_ode_eval:OutsideInterval');
bad = rod;
bad.qtilde_xi0(:) = Inf;
must_error_id(@() sph_radial_ode_eval(bad, ...
    radial_geometry(rod.type, rod.xi0)), ...
    'sph_radial_ode_eval:NonfiniteOutput');
end

function test_radial_weighted_endpoints
c = 0.8;
eg = sph_eigen('prolate', 1, c, 3);
g = radial_geometry('prolate', 1);
[V, ok] = sph_radial(eg, 1, g, c);
assert(ok);
assert(all(finite_radial_values(V)));

eg = sph_eigen('oblate', 4, c, 6);
g = radial_geometry('oblate', 0);
[V, ok] = sph_radial(eg, 1, g, c);
assert(ok);
assert(all(finite_radial_values(V)));

eg = sph_eigen('prolate', 0, c, 2);
g = radial_geometry('prolate', [1, 1.3]);
V = sph_radial(eg, 1, g, c);
assert(isequal(V.Va, zeros(size(V.V0))));
assert(isequal(V.Vm, zeros(size(V.V0))));
end

function test_radial_oblate_path_agreement
c = 0.9;
eg = sph_eigen('oblate', 3, c, 5);
h = 1e-5;
g = radial_geometry('oblate', [1 - h, 1, 1 + h]);
[V, ok] = sph_radial(eg, 1, g, c);
assert(ok);
names = {'V0', 'V1', 'V2', 'Va', 'Vm'};
for k = 1:numel(names)
    values = V.(names{k});
    assert_relative(values(:, 2), ...
        (values(:, 1) + values(:, 3))/2, 3e-9);
end
end

function test_radial_oblate_scaled_derivatives
mu = 40;
c = 1;
xi = 1.1e-3;
eg = unit_eigen('oblate', mu, mu, 0, c);
g = radial_geometry('oblate', xi);
[V, ok] = sph_radial(eg, 1, g, c);
assert(ok);

[T, dT, d2T] = scaled_bessel_series(mu, 0, c*xi);
alpha = g.alpha;
F = alpha^(mu/2);
L = mu*xi/alpha;
Lprime = mu*(1 - xi^2)/alpha^2;
expected = struct();
expected.V0 = F*c^mu*T;
expected.V1 = sqrt(alpha)*F*(c^(mu + 1)*dT + L*c^mu*T);
expected.V2 = alpha^(3/2)*F*(c^(mu + 2)*d2T + ...
    2*L*c^(mu + 1)*dT + (L^2 + Lprime)*c^mu*T);
expected.Va = expected.V0/sqrt(alpha);
expected.Vm = mu*expected.Va;
for name = {'V0', 'V1', 'V2', 'Va', 'Vm'}
    assert_relative(V.(name{1}), expected.(name{1}), 1e-10);
end
end

function test_radial_high_mu_combined_scale
mu = 100;
c = 2000;
xi = 1e-6;
eg = unit_eigen('oblate', mu, mu, 0, c);
g = radial_geometry('oblate', xi);
[V, ok] = sph_radial(eg, 1, g, c);
log_C = 0.5*log(pi) - (mu + 1)*log(2) - gammaln(mu + 1.5);
leading = exp(log_C + mu*log(c));
correction = 1 - (c*xi)^2/(2*(2*mu + 3));
assert(ok);
assert(isfinite(V.V0));
assert_relative(V.V0, g.alpha^(mu/2)*leading*correction, 3e-12);
end

function test_radial_high_mu_far_point
mu = 200;
c = 0.5;
eg = unit_eigen('prolate', mu, mu, 0, c);
[V, ok] = sph_radial(eg, 1, radial_geometry('prolate', 1e7), c);
assert(ok);
assert(all(finite_radial_values(V)));
end

function test_radial_high_c_regular_ivp
for c = [27.494593, 50, 200]
    eg = sph_eigen('prolate', 0, c, 2);
    xi = 1.03 + (0:8)/(100*c);
    g = radial_geometry('prolate', xi);
    [V, ok] = sph_radial(eg, 1, g, c);
    assert(ok && all(finite_radial_values(V)));
    assert_radial_ode(eg, V, g, c, 'prolate', 2e-11);

    h = xi(2) - xi(1);
    finite_difference = (V.V0(:, 1) - 8*V.V0(:, 2) + ...
        8*V.V0(:, 4) - V.V0(:, 5))/(12*h);
    returned = V.V1(:, 3)/sqrt(g.alpha(3));
    assert_relative(finite_difference, returned, 2e-5);

    separate = sph_radial(eg, 1, ...
        radial_geometry('prolate', xi(3)), c);
    assert_relative(separate.V0, V.V0(:, 3), 2e-11);
    [~, outgoing_ok] = sph_radial(eg, 3, ...
        radial_geometry('prolate', 3), c);
    assert(~outgoing_ok);
end
end

function test_radial_regular_ivp_scales
eg = ill_conditioned_eigen('prolate', 0, [0; 1], [0; 2], 0.5);
g = radial_geometry('prolate', 1 + (1:3)*1e-7);
[V, ok] = sph_radial(eg, 1, g, 0.5);
assert(ok && all(finite_radial_values(V)));

c = 20 + 1i;
eg = sph_eigen('prolate', 1, c, 2);
g = radial_geometry('prolate', [1, 1.02, 1.03]);
[V, ok] = sph_radial(eg, 1, g, c);
assert(ok && all(finite_radial_values(V)));
assert(all(V.V1(:, 1) == 1));
assert(all(V.V2(:, 1) == -1));
assert(all(V.Va(:, 1) == 1));
assert_radial_ode(eg, assign_points(V, 2:3), ...
    radial_geometry('prolate', [1.02, 1.03]), c, 'prolate', 2e-11);

eg = ill_conditioned_eigen('oblate', 0, [0; 1], [0; 2], 0.5);
[V, ok] = sph_radial(eg, 1, radial_geometry('oblate', [0, 0.2]), 0.5);
assert(ok && all(finite_radial_values(V)));
assert(isequal(V.V0(:, 1), [1; 0]));
assert(isequal(V.V1(:, 1), [0; 1]));
end

function test_radial_outgoing_tail_and_identity
c_tail = 50;
eg_tail = sph_eigen('prolate', 0, c_tail, 3);
g_near = radial_geometry('prolate', 1 + 1e-6);
[~, ok2] = sph_radial(eg_tail, 2, g_near, c_tail);
[~, ok3] = sph_radial(eg_tail, 3, g_near, c_tail);
assert(~ok2 && ~ok3);
must_error(@() sph_radial(eg_tail, 2, ...
    radial_geometry('prolate', 1), c_tail));

c = 0.5;
eg = sph_eigen('prolate', 0, c, 3);
g = radial_geometry('prolate', [2.5, 4]);
[V1, ok1] = sph_radial(eg, 1, g, c);
[V2, ok2] = sph_radial(eg, 2, g, c);
[V3, ok3] = sph_radial(eg, 3, g, c);
assert(ok1 && ok2 && ok3);
names = {'V0', 'V1', 'V2', 'Va', 'Vm'};
for k = 1:numel(names)
    assert_relative(V3.(names{k}), ...
        V1.(names{k}) + 1i*V2.(names{k}), 2e-12);
end
end

function test_radial_tail_parity_and_normalization
eg = unit_eigen('prolate', 0, 0, 10, 0.5);
eg.d(:, 1) = 0;
eg.d(1, 1) = 1;
eg.d(5, 1) = 1;
[~, ok] = sph_radial(eg, 2, radial_geometry('prolate', 4), 0.5);
assert(~ok);

eg = unit_eigen('prolate', 0, 0, 2, 0.5);
r = (0:2).';
log_weight = 0.5*(log(2*r + 1) + gammaln(r + 1) - gammaln(r + 1));
weight = exp(log_weight - max(log_weight));
eg.d(:, 1) = [1; 0; -weight(1)/weight(3)*(1 - eps/2)];
[~, ok] = sph_radial(eg, 2, radial_geometry('prolate', 2), 0.5);
assert(~ok);
end

function test_radial_outgoing_validation
eg = sph_eigen('prolate', 0, 0.5, 2);
g = radial_geometry('prolate', 2);
must_error(@() sph_radial(eg, 2, g, 0));
must_error(@() sph_radial(eg, 3, g, 0.5 + 0.01i));
must_error(@() sph_radial(eg, 1, g, 0.6));
end

function test_radial_validation_fails_closed
c = 20;
eg = sph_eigen('prolate', 0, c, 1);
g = radial_geometry('prolate', 1.03);

bad = eg;
bad.d(1) = NaN;
must_error_id(@() sph_radial(bad, 1, g, c), ...
    'sph_radial:InvalidEigenData');

bad = eg;
bad.lambda = bad.lambda(1);
must_error_id(@() sph_radial(bad, 1, g, c), ...
    'sph_radial:InvalidEigenData');

eg = sph_eigen('prolate', uint8(0), uint8(c), uint16(1));
[V, ok] = sph_radial(eg, 1, g, uint8(c));
assert(ok && all(finite_radial_values(V)));
end

function test_validation
must_error(@() sph_eigen('sphere', 0, 1, 2));
must_error(@() sph_eigen('prolate', -1, 1, 2));
must_error(@() sph_eigen('prolate', 1.5, 1, 2));
must_error(@() sph_eigen('prolate', 0, Inf, 2));
must_error(@() sph_eigen('prolate', 2, 1, 1));
must_error(@() sph_eigen('prolate', 0, 1, 401));
must_error_id(@() sph_eigen('prolate', 0, 1e6, 1), ...
    'sph_eigen:NoConvergence');
eg_integer = sph_eigen('prolate', uint8(0), uint8(200), uint16(2));
assert(isa(eg_integer.c, 'double') && eg_integer.c == 200);
must_error(@() sph_angular(struct(), struct('eta', 0, 'beta', 1)));
eg = sph_eigen('prolate', 0, 1, 1);
must_error(@() sph_angular(eg, struct('eta', [0; 1], 'beta', 1)));
must_error(@() sph_angular(eg, struct('eta', 2, 'beta', 0)));
end

function [T, dT, d2T] = scaled_bessel_series(mu, r, x)
T = complex(zeros(size(r)));
dT = T;
d2T = T;
for row = 1:numel(r)
    n = mu + r(row);
    C = inverse_odd_double_factorial(n);
    T(row) = C*x^r(row)*(1 - x^2/(2*(2*n + 3)) + ...
        x^4/(8*(2*n + 3)*(2*n + 5)));
    dT(row) = C*(r(row)*x^(max(r(row) - 1, 0)) - ...
        (r(row) + 2)*x^(r(row) + 1)/(2*(2*n + 3)) + ...
        (r(row) + 4)*x^(r(row) + 3)/ ...
        (8*(2*n + 3)*(2*n + 5)));
    if r(row) == 0
        dT(row) = C*(-x/(2*n + 3) + ...
            x^3/(2*(2*n + 3)*(2*n + 5)));
    end
    if r(row) >= 2
        leading = r(row)*(r(row) - 1)*x^(r(row) - 2);
    else
        leading = 0;
    end
    d2T(row) = C*(leading - ...
        (r(row) + 2)*(r(row) + 1)*x^r(row)/(2*(2*n + 3)) + ...
        (r(row) + 4)*(r(row) + 3)*x^(r(row) + 2)/ ...
        (8*(2*n + 3)*(2*n + 5)));
end
end

function value = inverse_odd_double_factorial(n)
value = exp(0.5*log(pi) - (n + 1)*log(2) - gammaln(n + 1.5));
end

function assert_radial_ode(eg, V, g, c, type, tolerance)
xi = g.xi(:).';
alpha = g.alpha(:).';
R = V.V0;
Rp = V.V1./sqrt(alpha);
Rpp = V.V2./alpha.^(3/2);
if strcmp(type, 'prolate')
    coefficient = eg.lambda - c^2*xi.^2 + eg.mu^2./alpha;
else
    coefficient = eg.lambda - c^2*xi.^2 - eg.mu^2./alpha;
end
terms = cat(3, alpha.*Rpp, 2*xi.*Rp, -coefficient.*R);
residual = sum(terms, 3);
scale = sum(abs(terms), 3);
assert(max(abs(residual(:))./max(scale(:), realmin)) < tolerance);
end

function residual = outgoing_fd_residual(rod, eg)
x = max(rod.xi0 + 0.4, 1.4);
h = min(5e-4, 0.01 / rod.c);
xi = x + h * (-2:2);
V = sph_radial_ode_eval(rod, radial_geometry(rod.type, xi));
R = V.V0;
Rp = (R(:, 1) - 8 * R(:, 2) + 8 * R(:, 4) - R(:, 5)) / ...
    (12 * h);
Rpp = (-R(:, 1) + 16 * R(:, 2) - 30 * R(:, 3) + ...
    16 * R(:, 4) - R(:, 5)) / (12 * h^2);
sigma = 1;
if strcmp(rod.type, 'oblate')
    sigma = -1;
end
alpha = x^2 - sigma;
terms = [alpha * Rpp, 2 * x * Rp, ...
    (rod.c^2 * x^2 - eg.lambda - sigma * eg.mu^2 / alpha) .* R(:, 3)];
residual = max(abs(sum(terms, 2)) ./ ...
    max(sum(abs(terms), 2), realmin));
end

function values = finite_radial_values(V)
values = [isfinite(V.V0(:)); isfinite(V.V1(:)); ...
    isfinite(V.V2(:)); isfinite(V.Va(:)); isfinite(V.Vm(:))];
end

function eg = unit_eigen(type, mu, l, rmax, c)
l = l(:);
d = zeros(rmax + 1, numel(l));
for column = 1:numel(l)
    d(l(column) - mu + 1, column) = 1;
end
eg = struct('type', type, 'mu', mu, 'c', c, 'L', max(l), ...
    'l', l, 'lambda', l.*(l + 1), 'd', d, ...
    'r', (0:rmax).', 'rmax', rmax);
end


function eg = ill_conditioned_eigen(type, mu, l, lambda, c)
rmax = 3;
r = (0:rmax).';
log_weight = 0.5*(log((2*mu + 2*r + 1)/(2*mu + 1)) + ...
    gammaln(r + 2*mu + 1) - gammaln(r + 1) - gammaln(2*mu + 1));
weight = exp(log_weight - max(log_weight));
d = zeros(rmax + 1, numel(l));
d(1, 1) = 1;
d(3, 1) = -weight(1)/weight(3)*(1 - eps/2);
d(2, 2) = 1;
d(4, 2) = -weight(2)/weight(4)*(1 - eps/2);
eg = struct('type', type, 'mu', mu, 'c', c, 'L', max(l), ...
    'l', l, 'lambda', lambda, 'd', d, 'r', r, 'rmax', rmax);
end

function selected = assign_points(V, points)
names = fieldnames(V);
for k = 1:numel(names)
    selected.(names{k}) = V.(names{k})(:, points);
end
end

function g = radial_geometry(type, xi)
g.xi = xi(:);
if strcmp(type, 'prolate')
    g.alpha = g.xi.^2 - 1;
else
    g.alpha = g.xi.^2 + 1;
end
end

function assert_relative(actual, expected, tolerance)
scale = max([norm(actual(:)), norm(expected(:)), realmin]);
assert(norm(actual(:) - expected(:))/scale <= tolerance);
end

function A = angular_matrix(type, mu, c, r)
ell = mu + r;
c2 = c^2;
if strcmp(type, 'oblate')
    c2 = -c2;
end
diagonal = ell .* (ell + 1) + ...
    (2 * ell .* (ell + 1) - 2 * mu^2 - 1) * c2 ./ ...
    ((2 * ell - 1) .* (2 * ell + 3));
rc = r(1:end-2);
off_diagonal = c2 * sqrt((2 * mu + rc + 2) .* (2 * mu + rc + 1) .* ...
    (rc + 2) .* (rc + 1)) ./ ...
    ((2 * mu + 2 * rc + 3) .* ...
    sqrt((2 * mu + 2 * rc + 1) .* (2 * mu + 2 * rc + 5)));
A = diag(diagonal) + diag(off_diagonal, 2) + diag(off_diagonal, -2);
end

function must_error(f)
did_error = false;
try
    f();
catch
    did_error = true;
end
assert(did_error);
end

function must_error_id(f, identifier)
try
    f();
catch exception
    assert(strcmp(exception.identifier, identifier));
    return
end
error('test_special_functions:ExpectedError', ...
    'Expected error %s was not raised.', identifier);
end
