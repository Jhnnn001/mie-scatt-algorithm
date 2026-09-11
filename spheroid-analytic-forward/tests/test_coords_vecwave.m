function test_coords_vecwave
% Coordinate-degeneracy and vector spheroidal wavefunction checks.

test_shapes_phase_and_m0;
test_generic_maxwell_identities(1);
test_generic_maxwell_identities(3);
test_generic_outgoing_ode_maxwell_identities;
test_degenerate_maxwell_identities;
test_regular_coordinate_sets;
test_solved_field_singularities;
test_validation;

fprintf('test_coords_vecwave: PASS\n');
end

function test_solved_field_singularities
semifocal = 1;
lambda = 4*pi;
n_m = 1;
n_p = 1.1;
theta_inc = 0;
phi_inc = 0;
pol = [1, 0.25i];
opts = struct('tol', 1e-4, 'Lmax', 60, 'Mmax', 40, ...
    'on_fail', 'error');
delta = 1e-8;
axis_offsets = [1e-6, 1e-7];
axis_z = [0.3; 2];

for type_cell = {'prolate', 'oblate'}
    type = type_cell{1};
    if strcmp(type, 'prolate')
        a = 2/sqrt(3);
        b = 1/sqrt(3);
    else
        a = 1/sqrt(3);
        b = 2/sqrt(3);
    end
    sol = spheroid_solve(lambda, n_p, n_m, a, b, theta_inc, ...
        phi_inc, pol, opts);
    assert(abs(sol.semifocal - semifocal) < 16*eps);
    [~, ~, ~, E0] = incident_plane_wave(sol.km, n_m, theta_inc, ...
        phi_inc, pol, 0, 0, 0);
    E0_scale = norm(E0);

    exact_axis = spheroid_eval(sol, zeros(size(axis_z)), ...
        zeros(size(axis_z)), semifocal*axis_z, 'total', 'auto');
    for rho = axis_offsets
        offset_axis = spheroid_eval(sol, semifocal*rho*ones(size(axis_z)), ...
            zeros(size(axis_z)), semifocal*axis_z, 'total', 'auto');
        assert(max(row_norm(offset_axis - exact_axis)) <= ...
            3*sol.c_ext*rho*E0_scale);
    end

    if strcmp(type, 'prolate')
        singular_points = semifocal*[0, 0, 1; delta/2, 0, 1];
    else
        singular_points = semifocal*[1, 0, 0; 1 + delta/2, 0, 0];
    end
    singular_field = spheroid_eval(sol, singular_points(:, 1), ...
        singular_points(:, 2), singular_points(:, 3), 'total', 'auto');
    assert(norm(singular_field(1, :) - singular_field(2, :)) <= ...
        1e-7*E0_scale);

    if strcmp(type, 'oblate')
        disk_points = semifocal*[0.45, 0, 1e-9; 0.45, 0, -1e-9];
        disk_field = spheroid_eval(sol, disk_points(:, 1), ...
            disk_points(:, 2), disk_points(:, 3), 'total', 'auto');
        assert(norm(disk_field(1, :) - disk_field(2, :)) <= ...
            1e-8*E0_scale);
    end
end
end

function value = row_norm(array)
value = sqrt(sum(abs(array).^2, 2));
end

function test_shapes_phase_and_m0
g = sph_coords('prolate', [0.8; 1.2], [0.3; -0.4], [0.5; 0.7]);
eg = sph_eigen('prolate', 0, 2.3, 2);
W = sph_angular(eg, g);
[V, ok] = sph_radial(eg, 1, g, 2.3);
assert(ok);
[M, N] = sph_vecwave(W, V, g, 0, 2.3);
assert_component_shapes(M, [3, 2]);
assert_component_shapes(N, [3, 2]);
assert(all(M.eta(:) == 0) && all(M.xi(:) == 0));
assert(all(N.phi(:) == 0));

g_shifted = g;
g_shifted.phi = g.phi + [0.4; -1.1];
[M_shifted, N_shifted] = sph_vecwave(W, V, g_shifted, 0, 2.3);
assert(isequaln(M, M_shifted));
assert(isequaln(N, N_shifted));

eg = sph_eigen('prolate', 2, 2.3, 4);
W = sph_angular(eg, g);
[V, ok] = sph_radial(eg, 1, g, 2.3);
assert(ok);
[M, N] = sph_vecwave(W, V, g, -2, 2.3);
g_shifted.phi = g.phi + [0.4; -1.1];
[M_shifted, N_shifted] = sph_vecwave(W, V, g_shifted, -2, 2.3);
assert(isequaln(M, M_shifted));
assert(isequaln(N, N_shifted));
end

function test_generic_maxwell_identities(kind)
c = 2.3;
for type_cell = {'prolate', 'oblate'}
    type = type_cell{1};
    if kind == 1
        if strcmp(type, 'prolate')
            point = [1.7, 0.8, 0.4];
        else
            point = [1.4, 0.6, 0.7];
        end
    else
        point = [4.0, 1.0, 0.8];
        assert(sph_coords(type, point(1), point(2), point(3)).xi > 2);
    end
    for m = [-3, 0, 2, 7]
        eg = sph_eigen(type, abs(m), c, abs(m) + 2);
        mode = 3;
        assert_cartesian_identities(type, eg, m, c, kind, point, mode, true);
    end
end
end

function test_generic_outgoing_ode_maxwell_identities
c = 2.3;
for type_cell = {'prolate', 'oblate'}
    type = type_cell{1};
    if strcmp(type, 'prolate')
        point = [4.0, 1.0, 0.8];
        xi0 = 1.2;
    else
        point = [4.0, 1.0, 0.8];
        xi0 = 0.5;
    end
    for m = [-3, 0, 2, 7]
        eg = sph_eigen(type, abs(m), c, abs(m) + 2);
        [rod, chk] = sph_radial_ode(eg, c, xi0, 1.5);
        assert(chk.ok);
        assert_cartesian_identities(type, eg, m, c, 3, point, 3, true, rod);
    end
end
end

function test_degenerate_maxwell_identities
c = 1.7;
cases = { ...
    'prolate', 1, [0, 0, 2]; ...
    'oblate', 1, [0, 0, 0.5]; ...
    'prolate', 1, [0, 0, 0.3]; ...
    'oblate', 1, [0.45, 0.2, 1e-9]; ...
    'oblate', 1, [0.45, 0.2, -1e-9]; ...
    'prolate', 3, [0, 0, 4]; ...
    'oblate', 3, [0, 0, 4]};
for row = 1:size(cases, 1)
    type = cases{row, 1};
    kind = cases{row, 2};
    point = cases{row, 3};
    m = 1;
    eg = sph_eigen(type, m, c, m + 2);
    assert_cartesian_identities(type, eg, m, c, kind, point, 3, false);
end
end

function assert_cartesian_identities(type, eg, m, c, kind, point, mode, check_scalar, rod)
if nargin < 9
    rod = [];
end
[M0, N0] = one_mode_cartesian(type, eg, m, c, kind, point, mode, rod);
h = 1e-5;
derivative_M = complex(zeros(3));
derivative_N = complex(zeros(3));
gradient_Pi = complex(zeros(1, 3));
for coordinate = 1:3
    step = zeros(1, 3);
    step(coordinate) = h;
    [Mplus, Nplus] = one_mode_cartesian( ...
        type, eg, m, c, kind, point + step, mode, rod);
    [Mminus, Nminus] = one_mode_cartesian( ...
        type, eg, m, c, kind, point - step, mode, rod);
    derivative_M(:, coordinate) = (Mplus - Mminus).'/(2*h);
    derivative_N(:, coordinate) = (Nplus - Nminus).'/(2*h);
    if check_scalar
        Pi_plus = scalar_wave(type, eg, m, c, kind, point + step, mode, rod);
        Pi_minus = scalar_wave(type, eg, m, c, kind, point - step, mode, rod);
        gradient_Pi(coordinate) = (Pi_plus - Pi_minus)/(2*h);
    end
end
if check_scalar
    assert(relative_error(M0, cross(gradient_Pi, point)) < 1e-6);
end
curl_M = matrix_curl(derivative_M);
curl_N = matrix_curl(derivative_N);
assert(relative_error(curl_M, c*N0) < 1e-6);
assert(relative_error(curl_N, c*M0) < 1e-6);
assert(abs(trace(derivative_M))/max(norm(derivative_M, 'fro'), realmin) < 1e-6);
assert(abs(trace(derivative_N))/max(norm(derivative_N, 'fro'), realmin) < 1e-6);
end

function test_regular_coordinate_sets
% Mode-level T2 regularity complements the solved-field T9 checks above.
c = 1.7;
for type_cell = {'prolate', 'oblate'}
    type = type_cell{1};
    for m = [-3, 0, 2]
        eg = sph_eigen(type, abs(m), c, abs(m) + 2);
        mode = 3;
        if strcmp(type, 'prolate')
            sets = {[0, 0, 2], [1e-7, 0, 2], ...
                [0, 0, 0.3], [1e-7, 0, 0.3]};
        else
            sets = {[0, 0, 0.5], [1e-7, 0, 0.5], ...
                [0.45, 0.2, 0], [0.45, 0.2, 1e-9], ...
                [0.45, 0.2, -1e-9]};
        end
        values = cell(size(sets));
        for k = 1:numel(sets)
            [M, N] = one_mode_cartesian(type, eg, m, c, 1, sets{k}, mode);
            values{k} = [M, N];
            assert(all(isfinite(real(values{k}))) && all(isfinite(imag(values{k}))));
        end
        if strcmp(type, 'prolate')
            assert_offset_limit(type, eg, m, c, mode, ...
                [0, 0, 2], values{1}, values{2});
            exact = values{3};
            near1 = values{4};
            [M2, N2] = one_mode_cartesian( ...
                type, eg, m, c, 1, [1e-6, 0, 0.3], mode);
            near2 = [M2, N2];
            assert(norm(exact - near1) <= 0.2*norm(exact - near2) + 1e-10);
        else
            assert_offset_limit(type, eg, m, c, mode, ...
                [0, 0, 0.5], values{1}, values{2});
            scale = max([norm(values{3}), norm(values{4}), norm(values{5}), 1]);
            assert(norm(values{3} - values{4})/scale < 1e-8);
            assert(norm(values{4} - values{5})/scale < 1e-8);
        end
    end
end

% Compare the effective delta and delta/2 singular-set nudges.
delta = 1e-8;
for item = {{'prolate', [0, 0, 1], [delta/2, 0, 1]}, ...
        {'oblate', [1, 0, 0], [1 + delta/2, 0, 0]}}
    type = item{1}{1};
    point = item{1}{2};
    half_delta_point = item{1}{3};
    eg = sph_eigen(type, 2, c, 4);
    [M, N] = one_mode_cartesian(type, eg, 2, c, 1, point, 3);
    [Mh, Nh] = one_mode_cartesian(type, eg, 2, c, 1, half_delta_point, 3);
    values = [M, N];
    half_delta_values = [Mh, Nh];
    assert(all(isfinite(real([values, half_delta_values]))) && ...
        all(isfinite(imag([values, half_delta_values]))));
    scale = max([norm(values), norm(half_delta_values), 1]);
    assert(norm(values - half_delta_values)/scale <= 1e-7);
end
end

function assert_offset_limit(type, eg, m, c, mode, point, exact, near1)
point(1) = 1e-6;
[M2, N2] = one_mode_cartesian(type, eg, m, c, 1, point, mode);
near2 = [M2, N2];
assert(norm(exact - near1) <= 0.2*norm(exact - near2) + 1e-10);
end

function [Mcart, Ncart] = one_mode_cartesian(type, eg, m, c, kind, point, mode, rod)
if nargin < 8
    rod = [];
end
g = sph_coords(type, point(1), point(2), point(3));
W = sph_angular(eg, g);
if isempty(rod)
    [V, ok] = sph_radial(eg, kind, g, c);
    assert(ok);
else
    V = sph_radial_ode_eval(rod, g);
end
[M, N] = sph_vecwave(W, V, g, m, c);
phase = exp(1i*m*g.phi);
Mcart = phase*(M.eta(mode)*g.e_eta + M.xi(mode)*g.e_xi + ...
    M.phi(mode)*g.e_phi);
Ncart = phase*(N.eta(mode)*g.e_eta + N.xi(mode)*g.e_xi + ...
    N.phi(mode)*g.e_phi);
end

function Pi = scalar_wave(type, eg, m, c, kind, point, mode, rod)
if nargin < 8
    rod = [];
end
g = sph_coords(type, point(1), point(2), point(3));
W = sph_angular(eg, g);
if isempty(rod)
    [V, ok] = sph_radial(eg, kind, g, c);
    assert(ok);
else
    V = sph_radial_ode_eval(rod, g);
end
Pi = W.W0(mode)*V.V0(mode)*exp(1i*m*g.phi);
end

function curl_value = matrix_curl(D)
curl_value = [D(3, 2) - D(2, 3), ...
    D(1, 3) - D(3, 1), D(2, 1) - D(1, 2)];
end

function value = relative_error(actual, expected)
value = norm(actual - expected)/max([norm(actual), norm(expected), realmin]);
end

function assert_component_shapes(value, expected)
assert(isequal(sort(fieldnames(value)), {'eta'; 'phi'; 'xi'}));
for name = {'eta', 'xi', 'phi'}
    assert(isequal(size(value.(name{1})), expected));
end
end

function test_validation
g = sph_coords('prolate', 1, 0.2, 0.3);
eg = sph_eigen('prolate', 1, 2, 2);
W = sph_angular(eg, g);
[V, ok] = sph_radial(eg, 1, g, 2);
assert(ok);
must_error_id(@() sph_vecwave(struct(), V, g, 1, 2), ...
    'sph_vecwave:InvalidAngularData');
must_error_id(@() sph_vecwave(W, struct(), g, 1, 2), ...
    'sph_vecwave:InvalidRadialData');
must_error_id(@() sph_vecwave(W, V, rmfield(g, 'D'), 1, 2), ...
    'sph_vecwave:InvalidCoordinates');
must_error_id(@() sph_vecwave(W, V, g, 0.5, 2), ...
    'sph_vecwave:InvalidAzimuthalIndex');
must_error_id(@() sph_vecwave(W, V, g, 1, 0), ...
    'sph_vecwave:InvalidParameter');
bad = W;
bad.W0 = [bad.W0, bad.W0];
must_error_id(@() sph_vecwave(bad, V, g, 1, 2), ...
    'sph_vecwave:SizeMismatch');
bad = g;
bad.B_prime = bad.B_prime + 0.25;
must_error_id(@() sph_vecwave(W, V, bad, 1, 2), ...
    'sph_vecwave:InvalidCoordinates');
end

function must_error_id(f, identifier)
try
    f();
catch exception
    assert(strcmp(exception.identifier, identifier));
    return
end
error('test_coords_vecwave:ExpectedError', ...
    'Expected error %s was not raised.', identifier);
end
