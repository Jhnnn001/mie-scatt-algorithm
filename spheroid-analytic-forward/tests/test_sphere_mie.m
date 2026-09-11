function test_sphere_mie
% Sphere Mie coefficients, fields, covariance, and interface conditions.

test_outside_envelope_warning;
test_efficiencies;
test_field_origin_and_axes;
test_no_contrast_field;
test_incident_without_solve;
test_rotation_and_modes;
test_boundary_conditions;
test_near_sphere_dispatch;
test_near_sphere;
fprintf('test_sphere_mie: PASS\n');
end

function test_incident_without_solve
lambda = 2*pi;
n_p = 1.3;
n_m = 1;
a = 1;
b = 0.75;
theta_inc = 0.4;
phi_inc = -0.3;
pol = [1, 0.2i];
X = [0.1; 1.4];
Y = [-0.2; 0.3];
Z = [0.25; -0.8];
opts = struct('tol', 1e-12, 'field', 'incident', 'side', 'auto', ...
    'Lmax', 1, 'Mmax', 1);
[Ex, Ey, Ez] = spheroid_field(lambda, n_p, n_m, a, b, ...
    theta_inc, phi_inc, pol, X, Y, Z, opts);
[Eincident, ~] = incident_plane_wave(2*pi*n_m/lambda, n_m, ...
    theta_inc, phi_inc, pol, X, Y, Z);
assert(isequal([Ex(:), Ey(:), Ez(:)], Eincident));

opts.side = 'out';
must_error_id(@() spheroid_field(lambda, n_p, n_m, a, b, ...
    theta_inc, phi_inc, pol, 0, 0, 2*a, opts), ...
    'spheroid_field:InvalidForcedSide');
end

function test_near_sphere_dispatch
b = 1;
a = b*(1 + 5e-7);
lambda = 100*pi;
n_p = 1.02;
n_m = 1;
theta_inc = 0;
phi_inc = 0;
pol = [1, 0];
opts = struct('tol', 3e-3, 'field', 'total', 'side', 'auto', ...
    'Lmax', 30, 'Mmax', 20, 'on_fail', 'error', 'verbose', false);
[Ex, Ey, Ez, info] = spheroid_field(lambda, n_p, n_m, a, b, ...
    theta_inc, phi_inc, pol, b, 0, 0, opts);
assert(strcmp(info.route, 'spheroid'));
[Eout, ~] = spheroid_eval(info.sol, b, 0, 0, 'total', 'out');
assert_relative([Ex, Ey, Ez], Eout, 5e-13);
empty_mu = cellfun(@isempty, info.sol.mu_data);
assert(isequal(find(~empty_mu), 2));
assert(numel(info.rmax_by_mu.external) == info.M + 1);
assert(all(isnan(info.rmax_by_mu.external(empty_mu))));
assert(info.rmax_by_mu.external(2) == info.sol.mu_data{2}.eg_ext.rmax);
assert(isfinite(info.rcond_min) && info.rcond_min == info.rcond_min_active);
skipped = abs([info.sol.modes.m]) ~= 1;
assert(all([info.sol.modes(skipped).analytically_zeroed]));
assert(all(cellfun(@(value) all(value == 0), ...
    {info.sol.modes(skipped).a, info.sol.modes(skipped).b, ...
    info.sol.modes(skipped).c, info.sol.modes(skipped).d})));
bad = info.sol;
bad.mu_data{2} = [];
must_error_id(@() spheroid_eval(bad, b, 0, 0, 'total', 'out'), ...
    'spheroid_eval:InvalidSolution');
bad = info.sol;
bad.mu_data = struct();
must_error_id(@() spheroid_eval(bad, b, 0, 0, 'total', 'out'), ...
    'spheroid_eval:InvalidSolution');
bad = info.sol;
bad.mu_data{2} = rmfield(bad.mu_data{2}, 'eg_ext');
must_error_id(@() spheroid_eval(bad, b, 0, 0, 'total', 'out'), ...
    'spheroid_eval:InvalidSolution');
bad = info.sol;
bad.mu_data{2}.eg_ext.c = bad.mu_data{2}.eg_ext.c + 1;
must_error_id(@() spheroid_eval(bad, b, 0, 0, 'total', 'out'), ...
    'spheroid_eval:InvalidSolution');
bad = info.sol;
index = find(skipped, 1);
bad.modes(index).a = [];
must_error_id(@() spheroid_eval(bad, b, 0, 0, 'total', 'out'), ...
    'spheroid_eval:InvalidSolution');

back = spheroid_solve(lambda, n_p, n_m, a, b, pi, phi_inc, pol, opts);
assert(isequal(find(~cellfun(@isempty, back.mu_data)), 2));
end

function test_near_sphere
n_m = 1;
n_p = 1.2;
b = 1;
lambda = 2*pi/10;
axis_values = linspace(-1.5*b, 1.5*b, 41);
[X, Z] = meshgrid(axis_values, axis_values);
Y = 0.3*b*ones(size(X));
cases = struct( ...
    'theta', {0, 40*pi/180}, ...
    'phi', {0, 70*pi/180}, ...
    'pol', {[1, 0], [1, 0.5i]});
opts = struct('tol', 1e-4, 'field', 'total', 'side', 'auto', ...
    'Lmax', 80, 'Mmax', 60, 'on_fail', 'error', 'verbose', false);
for aspect = [1 + 1e-4, 1 - 1e-4]
    a = aspect*b;
    for k = 1:numel(cases)
        illumination = cases(k);
        [Ex, Ey, Ez, info] = spheroid_field(lambda, n_p, n_m, a, b, ...
            illumination.theta, illumination.phi, illumination.pol, ...
            X, Y, Z, opts);
        spheroid_E = [Ex(:), Ey(:), Ez(:)];
        [sphere_E, ~] = mie_field(lambda, n_p, n_m, b, ...
            illumination.theta, illumination.phi, illumination.pol, ...
            X, Y, Z, 'total', 'auto');
        [~, ~, ~, E0] = incident_plane_wave(2*pi*n_m/lambda, n_m, ...
            illumination.theta, illumination.phi, illumination.pol, ...
            0, 0, 0);
        error_value = max(sqrt(sum(abs(spheroid_E - sphere_E).^2, 2)))/ ...
            norm(E0);
        assert(error_value < 1e-2);
        assert(strcmp(info.route, 'spheroid') && info.validated);
        assert_validation_details(info);
        if illumination.theta ~= 0
            modes = info.sol.modes;
            assert(~any([modes.roundoff_zeroed] & [modes.rhs_norm] > 0));
        end
        active = ~[info.sol.modes.analytically_zeroed];
        assert(all(isfinite([info.sol.modes(active).numerical_rank])));
        assert(all(ismember({info.sol.modes.solver_method}, ...
            {'analytic_zero', 'lsqminnorm', 'backslash'})));
        if aspect > 1 && k == 1
            assert_far_scattered_scaling(info.sol, b);
        end
        Xboundary = [0; b/sqrt(2); 0.6*b];
        Yboundary = [0; b/sqrt(2); 0];
        Zboundary = [a; 0; 0.8*a];
        [Eauto, Hauto] = spheroid_eval(info.sol, Xboundary, ...
            Yboundary, Zboundary, ...
            'total', 'auto');
        [Eout, Hout] = spheroid_eval(info.sol, Xboundary, ...
            Yboundary, Zboundary, ...
            'total', 'out');
        assert_relative(Eauto, Eout, 2e-12);
        assert_relative(Hauto, Hout, 2e-12);
    end
end

[Ex, Ey, Ez, info] = spheroid_field(lambda, n_p, n_m, b, b, ...
    cases(2).theta, cases(2).phi, cases(2).pol, X, Y, Z, opts);
[sphere_E, ~] = mie_field(lambda, n_p, n_m, b, cases(2).theta, ...
    cases(2).phi, cases(2).pol, X, Y, Z, 'total', 'auto');
assert(strcmp(info.route, 'mie'));
assert(isequal([Ex(:), Ey(:), Ez(:)], sphere_E));

threshold_a = b + 32*eps(b);
[~, ~, ~, info] = spheroid_field(lambda, n_m, n_m, threshold_a, b, ...
    cases(2).theta, cases(2).phi, cases(2).pol, 0, 0, 0, opts);
assert(strcmp(info.route, 'mie'));

tight_opts = opts;
tight_opts.tol = 1e-12;
[~, ~, ~, info] = spheroid_field(lambda, n_p, n_m, b, b, ...
    cases(2).theta, cases(2).phi, cases(2).pol, 0, 0, 0, tight_opts);
assert(~info.validated && ~info.checks.mie_tolerance.passed);
end

function assert_far_scattered_scaling(sol, scale)
radius = scale*[1e6; 1e18];
E = spheroid_eval(sol, zeros(2, 1), zeros(2, 1), radius, ...
    'scattered', 'auto');
scaled_amplitude = radius.*sqrt(sum(abs(E(:, 1:2)).^2, 2));
assert(all(isfinite(scaled_amplitude)) && all(scaled_amplitude > 0));
assert(abs(scaled_amplitude(2)/scaled_amplitude(1) - 1) < 1e-4);
end

function test_outside_envelope_warning
identifier = 'spheroid_field:OutsideValidatedEnvelope';
state = warning('query', identifier);
restore_warning = onCleanup(@() warning(state.state, identifier));
warning('on', identifier);
solver_identifier = 'spheroid_solve:OutsideValidatedEnvelope';
solver_state = warning('query', solver_identifier);
restore_solver_warning = onCleanup( ...
    @() warning(solver_state.state, solver_identifier));
warning('on', solver_identifier);
lastwarn('');
[~, ~, ~, info] = spheroid_field(2*pi/300, 1, 1, 1, 0.5, ...
    0.3, 0.4, [1, 0], 0, 0, 0, struct());
[~, actual_identifier] = lastwarn;
assert(strcmp(actual_identifier, identifier));
assert(~info.validated);
lastwarn('');
[~, ~, ~, info] = spheroid_field(200*pi, 1.01, 1, 1, 0.5, ...
    0.3, 0.4, [1, 0], 0, 0, 0, struct());
[~, actual_identifier] = lastwarn;
assert(strcmp(actual_identifier, solver_identifier));
assert(~info.validated);
end

function assert_validation_details(info)
names = fieldnames(info.checks);
for k = 1:numel(names)
    check = info.checks.(names{k});
    assert(check.computed && check.passed);
end
surface = info.checks.surface_residual.details;
assert(surface.refined.maximum <= ...
    1.1*max(surface.coarse.maximum, 64*eps));
for k = 1:numel(info.special_function_diagnostics)
    radial = info.special_function_diagnostics(k).radial;
    assert(radial.ok);
    assert(radial.ode_residual <= radial.backward_error_limit);
    assert(radial.asymptotic_residual <= radial.backward_error_limit);
    assert(radial.ode_residual_floor <= radial.backward_error_limit);
    assert(radial.asymptotic_residual_floor <= radial.backward_error_limit);
end
end

function test_no_contrast_field
lambda = 0.8;
n_m = 1.2;
a = 0.5;
theta_inc = 0.4;
phi_inc = 0.7;
pol = [1, 0.2i];
points = a*[0.2, 0.3, 0.4; -0.6, 0.1, 0.2];
[E, H] = mie_field(lambda, n_m, n_m, a, theta_inc, phi_inc, pol, ...
    points(:, 1), points(:, 2), points(:, 3), 'total');
[Ei, Hi] = incident_plane_wave(2*pi*n_m/lambda, n_m, theta_inc, ...
    phi_inc, pol, points(:, 1), points(:, 2), points(:, 3));
assert_relative(E, Ei, 2e-12);
assert_relative(H, Hi, 2e-12);
end

function test_efficiencies
% Wiscombe MVTstNew case 14 uses m=1.5-1i with the opposite absorption
% convention.  For exp(-i*omega*t), use m=1.5+1i.  Source:
% https://raw.githubusercontent.com/cfinch/Mie_scattering/master/Wiscombe/MVTstNew.f
[qext, qsca, qback] = mie_efficiencies(1, 1.5 + 1i);
assert_close(qext, 2.336321, 2e-6);
assert_close(qsca, 0.6634538, 2e-6);
assert_close(qback, 0.573002538306, 2e-6);

[qext, qsca, qback, coeff] = mie_efficiencies(3.2, 1);
assert(isequal(coeff.n, (1:numel(coeff.n)).'));
assert(isequal(coeff.a, zeros(size(coeff.a))));
assert(isequal(coeff.b, zeros(size(coeff.b))));
assert(isequal(coeff.c, ones(size(coeff.c))));
assert(isequal(coeff.d, ones(size(coeff.d))));
assert(qext == 0 && qsca == 0 && qback == 0);

x = 0.01;
m = 1.4;
[~, qsca] = mie_efficiencies(x, m);
rayleigh = (8/3)*x^4*abs((m^2 - 1)/(m^2 + 2))^2;
assert(abs(qsca - rayleigh)/rayleigh < 1e-4);

[qext, qsca, ~, coeff] = mie_efficiencies(4.7, 1.33);
w = 2*coeff.n + 1;
optical_residual = sum(w .* (real(coeff.a + coeff.b) - ...
    abs(coeff.a).^2 - abs(coeff.b).^2));
assert(abs(optical_residual) < 2e-11*sum(w));
assert(abs(qext - qsca) < 2e-12*max(1, qext));

[qext, qsca] = mie_efficiencies(8.3, 1.45 + 0.08i);
assert(qext - qsca >= -5e-13*max(1, qext));

[~, ~, ~, coeff] = mie_efficiencies(0.02, 1.4);
assert(imag(coeff.a(1)) < 0);
must_error_id(@() mie_efficiencies(1e-25, 1.5), ...
    'mie_efficiencies:NoConvergence');
end

function test_field_origin_and_axes
lambda = 0.83;
n_p = 1.45 + 0.03i;
n_m = 1.1;
a = 0.41;
theta_inc = 0.7;
phi_inc = -0.4;
pol = [0.8 + 0.2i, -0.3 + 0.5i];
x = 2*pi*n_m*a/lambda;
[~, ~, ~, coeff] = mie_efficiencies(x, n_p/n_m);
R = incidence_frame(theta_inc, phi_inc);

[E0, H0] = mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, 0, 0, 0, 'total');
expected_E0 = coeff.d(1)*([pol, 0]*R.');
expected_H0 = n_p*coeff.c(1)*([-pol(2), pol(1), 0]*R.');
assert_relative(E0, expected_E0, 5e-13);
assert_relative(H0, expected_H0, 5e-13);

z = [-1.6*a; 1.6*a];
eps_axis = 1e-8*a;
[E_axis, H_axis] = mie_field(lambda, n_p, n_m, a, 0, 0, ...
    [1, 0.4i], zeros(2, 1), zeros(2, 1), z, 'scattered');
[E_near, H_near] = mie_field(lambda, n_p, n_m, a, 0, 0, ...
    [1, 0.4i], eps_axis*ones(2, 1), zeros(2, 1), z, 'scattered');
assert(all(isfinite(E_axis(:))) && all(isfinite(H_axis(:))));
assert_relative(E_axis, E_near, 2e-7);
assert_relative(H_axis, H_near, 2e-7);

[E_axis, H_axis] = mie_field(lambda, n_p, n_m, a, 0, 0, ...
    [1, 0.4i], 0, 0, 0.35*a, 'total');
[E_near, H_near] = mie_field(lambda, n_p, n_m, a, 0, 0, ...
    [1, 0.4i], eps_axis, 0, 0.35*a, 'total');
assert(all(isfinite(E_axis(:))) && all(isfinite(H_axis(:))));
assert_relative(E_axis, E_near, 2e-7);
assert_relative(H_axis, H_near, 2e-7);
end

function test_rotation_and_modes
lambda = 0.74;
n_p = 1.38 + 0.02i;
n_m = 1.07;
a = 0.32;
theta_inc = 0.9;
phi_inc = 1.2;
pol = [0.7 - 0.1i, -0.2 + 0.6i];
R = incidence_frame(theta_inc, phi_inc);
points_local = a*[0.2, -0.1, 0.3; 1.2, 0.4, -0.7; -0.5, 0.8, 1.1];
points_world = points_local*R.';

[E_local, H_local] = mie_field(lambda, n_p, n_m, a, 0, 0, pol, ...
    points_local(:, 1), points_local(:, 2), points_local(:, 3), 'scattered');
[E_world, H_world] = mie_field(lambda, n_p, n_m, a, ...
    theta_inc, phi_inc, pol, points_world(:, 1), points_world(:, 2), ...
    points_world(:, 3), 'scattered');
assert_relative(E_world, E_local*R.', 2e-12);
assert_relative(H_world, H_local*R.', 2e-12);

X = [0; 1.3*a; a];
Y = [0.1*a; -0.2*a; 0];
Z = [0; 0.3*a; 0];
[Ei, Hi] = mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, X, Y, Z, 'incident');
[Ei_ref, Hi_ref] = incident_plane_wave(2*pi*n_m/lambda, n_m, ...
    theta_inc, phi_inc, pol, X, Y, Z);
[Et, Ht] = mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, X, Y, Z, 'total');
[Es, Hs] = mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, X, Y, Z, 'scattered');
assert(isequal(size(Et), [numel(X), 3]));
assert(isequal(size(Ht), [numel(X), 3]));
assert_relative(Ei, Ei_ref, 5e-15);
assert_relative(Hi, Hi_ref, 5e-15);
assert_relative(Et - Es, Ei, 2e-12);
assert_relative(Ht - Hs, Hi, 2e-12);

[E_auto, H_auto] = mie_field(lambda, n_p, n_m, a, theta_inc, ...
    phi_inc, pol, a, 0, 0, 'total');
[E_out, H_out] = mie_field(lambda, n_p, n_m, a, theta_inc, ...
    phi_inc, pol, a, 0, 0, 'total', 'out');
assert_relative(E_auto, E_out, 5e-14);
assert_relative(H_auto, H_out, 5e-14);

r_tie = a*(1 - 32*eps);
[E_auto, H_auto] = mie_field(lambda, n_p, n_m, a, theta_inc, ...
    phi_inc, pol, r_tie, 0, 0, 'total', 'auto');
[E_out, H_out] = mie_field(lambda, n_p, n_m, a, theta_inc, ...
    phi_inc, pol, r_tie, 0, 0, 'total', 'out');
assert(isequal(E_auto, E_out));
assert(isequal(H_auto, H_out));
must_error(@() mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, 0, 0, 0, 'total', 'out'));
must_error(@() mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, 2*a, 0, 0, 'total', 'in'));
end

function test_boundary_conditions
lambda = 0.91;
n_p = 1.52 + 0.04i;
n_m = 1.13;
a = 0.37;
theta_inc = 0.6;
phi_inc = -0.8;
pol = [1, 0.35i];
theta = [0.4; 1.1; 2.3];
phi = [-0.7; 0.2; 1.8];
er = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
points = a*er;
[Ein, Hin] = mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, points(:, 1), points(:, 2), points(:, 3), 'total', 'in');
[Eout, Hout] = mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, points(:, 1), points(:, 2), points(:, 3), 'total', 'out');

Ein_r = sum(Ein.*er, 2);
Eout_r = sum(Eout.*er, 2);
Hin_r = sum(Hin.*er, 2);
Hout_r = sum(Hout.*er, 2);
Ein_t = Ein - Ein_r.*er;
Eout_t = Eout - Eout_r.*er;
Hin_t = Hin - Hin_r.*er;
Hout_t = Hout - Hout_r.*er;
assert_relative(Ein_t, Eout_t, 3e-11);
assert_relative(Hin_t, Hout_t, 3e-11);
assert_relative(n_p^2*Ein_r, n_m^2*Eout_r, 3e-11);
assert_relative(Hin_r, Hout_r, 3e-11);
end

function R = incidence_frame(theta, phi)
e_te = [-sin(phi), cos(phi), 0];
e_tm = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)];
e1 = cos(phi)*e_tm - sin(phi)*e_te;
e2 = sin(phi)*e_tm + cos(phi)*e_te;
khat = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
R = [e1(:), e2(:), khat(:)];
end

function assert_close(actual, expected, tolerance)
assert(max(abs(actual(:) - expected(:))) <= ...
    tolerance*max(1, max(abs(expected(:)))));
end

function assert_relative(actual, expected, tolerance)
scale = max([norm(actual(:)), norm(expected(:)), 1e-14]);
assert(norm(actual(:) - expected(:))/scale <= tolerance);
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
did_error = false;
try
    f();
catch exception
    did_error = strcmp(exception.identifier, identifier);
end
assert(did_error);
end
