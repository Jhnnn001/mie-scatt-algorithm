function test_unit_scaling
% Shared length-unit invariance at extreme finite scales.

test_mie_shared_units;
test_spheroid_shared_units;
fprintf('test_unit_scaling: PASS\n');
end

function test_mie_shared_units
lambda = 2*pi;
n_p = 1.4 + 0.03i;
n_m = 1;
a = 1;
theta_inc = 0.6;
phi_inc = -0.4;
pol = [1, 0.2i];
points = [0.2, -0.3, 0.4; 1.8, 0.2, -0.4; -1.2, 0.7, 0.5];
[Eref, Href] = mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, points(:, 1), points(:, 2), points(:, 3), 'total', 'auto');
for factor = [1e-200, 1e200]
    scaled = factor*points;
    [E, H] = mie_field(factor*lambda, n_p, n_m, factor*a, ...
        theta_inc, phi_inc, pol, scaled(:, 1), scaled(:, 2), ...
        scaled(:, 3), 'total', 'auto');
    assert_relative(E, Eref, 2e-12);
    assert_relative(H, Href, 2e-12);
end
end

function test_spheroid_shared_units
lambda = 100*pi;
n_p = 1.01;
n_m = 1;
a = 1;
b = 0.75;
theta_inc = 0;
phi_inc = 0.3;
pol = [1, 0];
points = [0.1, 0.2, 0.15; 1.2, 0.2, -0.3];
opts = struct('tol', 3e-3, 'Lmax', 25, 'Mmax', 20, ...
    'on_fail', 'error', 'verbose', false);
factors = [1e-200, 1e200];
amplitudes = [1e-200, 1e200];
fields = cell(size(factors));
for k = 1:numel(factors)
    factor = factors(k);
    sol = spheroid_solve(factor*lambda, n_p, n_m, factor*a, ...
        factor*b, theta_inc, phi_inc, amplitudes(k)*pol, opts);
    assert(sol.info.validated);
    assert(abs(sol.xi0 - a/sqrt(a^2 - b^2)) < 32*eps);
    scaled = factor*points;
    [E, H] = spheroid_eval(sol, scaled(:, 1), scaled(:, 2), ...
        scaled(:, 3), 'total', 'auto');
    fields{k} = [E, H]/amplitudes(k);
end
assert_relative(fields{2}, fields{1}, 2e-10);
end

function assert_relative(actual, expected, tolerance)
scale = max([norm(actual(:)), norm(expected(:)), realmin]);
assert(norm(actual(:) - expected(:))/scale <= tolerance);
end
