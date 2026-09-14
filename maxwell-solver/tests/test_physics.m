function test_physics
%TEST_PHYSICS Physical invariants for the homogeneous-spheroid solver.
% T6 and T7 use the 2:1 prolate/oblate pair from T4.  The Mie far-field
% radius is 1e7*a because a sphere has zero semifocal length.

previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng));
rng(20260912, 'twister');

test_zero_contrast;
test_rayleigh;
solutions = test_energy;
test_boundary_traces(solutions);
fprintf('test_physics: PASS\n');
end

function test_zero_contrast
lambda = 0.8;
n_m = 1.2;
theta_inc = 35*pi/180;
phi_inc = 0.4;
pol = [0.8 + 0.1i, -0.3 + 0.2i];
k_m = 2*pi*n_m/lambda;
[~, ~, ~, E0] = incident_plane_wave(k_m, n_m, theta_inc, ...
    phi_inc, pol, 0, 0, 0);

for type = {'prolate', 'oblate'}
    [a, b] = axes_for_type(type{1}, 0.5);
    cloud = (2*rand(100, 3) - 1).*[1.6*b, 1.6*b, 1.6*a];
    opts = struct('tol', 1e-6, 'field', 'scattered');
    [Ex, Ey, Ez, info] = spheroid_field(lambda, n_m, n_m, a, b, ...
        theta_inc, phi_inc, pol, cloud(:, 1), cloud(:, 2), ...
        cloud(:, 3), opts);
    Escattered = [Ex, Ey, Ez];
    assert(max(row_norm(Escattered)) <= 1e-8*norm(E0));
    assert(strcmp(info.route, 'zero_contrast') && info.validated);

    direction = randn(100, 3);
    direction = direction./row_norm(direction);
    radius = 0.9*rand(100, 1).^(1/3);
    inside = radius.*direction.*[b, b, a];
    opts.field = 'total';
    [Ex, Ey, Ez] = spheroid_field(lambda, n_m, n_m, a, b, ...
        theta_inc, phi_inc, pol, inside(:, 1), inside(:, 2), ...
        inside(:, 3), opts);
    Etotal = [Ex, Ey, Ez];
    Eincident = incident_plane_wave(k_m, n_m, theta_inc, phi_inc, ...
        pol, inside(:, 1), inside(:, 2), inside(:, 3));
    assert_max_field(Etotal, Eincident, norm(E0), 1e-8);
end
end

function test_rayleigh
lambda = 2*pi;
n_m = 1;
n_p = 1.5;
x = 0.02;
opts = struct('tol', 1e-6, 'on_fail', 'error');
paired_coordinates = [0.28, 0.09, 0.13; ...
                     -0.17, 0.24, 0.11; ...
                      0.08, -0.19, 0.26];

for type = {'prolate', 'oblate'}
    [a, b] = axes_for_type(type{1}, x);
    [Lx, Lz] = depolarization_factors(type{1}, a, b);
    for orientation = 1:2
        if orientation == 1
            theta_inc = 0;
        else
            theta_inc = pi/2;
        end
        pol = [1, 0];
        sol = spheroid_solve(lambda, n_p, n_m, a, b, theta_inc, ...
            0, pol, opts);
        assert(sol.info.validated);
        if orientation == 1
            axial = abs([sol.modes.m]) ~= 1;
            assert(all([sol.modes(axial).analytically_zeroed]));
        elseif strcmp(type{1}, 'prolate')
            assert(any(strcmp({sol.modes.solver_method}, 'lsqminnorm')));
            active_rank = [sol.modes.numerical_rank];
            assert(all(isfinite(active_rank(...
                ~[sol.modes.analytically_zeroed]))));
        end
        [~, ~, ~, E0] = incident_plane_wave(sol.km, n_m, ...
            theta_inc, 0, pol, 0, 0, 0);
        epsilon_relative = (n_p/n_m)^2;
        expected = E0./(1 + [Lx, Lx, Lz]*(epsilon_relative - 1));

        [Ecenter, ~] = spheroid_eval(sol, 0, 0, 0, 'total', 'auto');
        assert_relative(Ecenter, expected, 1e-3);

        points = paired_coordinates.*[b, b, a];
        [Eplus, ~] = spheroid_eval(sol, points(:, 1), points(:, 2), ...
            points(:, 3), 'total', 'auto');
        [Eminus, ~] = spheroid_eval(sol, -points(:, 1), -points(:, 2), ...
            -points(:, 3), 'total', 'auto');
        % Inversion averaging removes the O(k*r) plane-wave phase gradient,
        % leaving the O(x^2) correction compared by the Rayleigh criterion.
        even_field = (Eplus + Eminus)/2;
        expected_points = repmat(expected, size(points, 1), 1);
        assert_max_field(even_field, expected_points, norm(E0), 1e-3);
    end
end
end

function solutions = test_energy
calibrate_far_field_with_mie;

lambda = 2*pi;
n_m = 1;
x = 30;
theta_inc = 30*pi/180;
phi_inc = 0;
pol = [0, 1];
opts = struct('tol', 1e-6, 'on_fail', 'error');
energy_tolerance = 1e-5;
quadrature_tolerance = energy_tolerance/4;
solutions = cell(2, 1);

for k = 1:2
    types = {'prolate', 'oblate'};
    type = types{k};
    [a, b] = axes_for_type(type, x);

    lossless = spheroid_solve(lambda, 1.5, n_m, a, b, ...
        theta_inc, phi_inc, pol, opts);
    assert(lossless.info.validated);
    [Csca, Cext] = far_cross_sections(lossless);
    assert_relative(Cext, Csca, 1e-6);
    fprintf('T7 %s lossless: Cext %.12g, Csca %.12g, rel %.3g\n', ...
        type, Cext, Csca, relative_gap(Cext, Csca));

    lossy = spheroid_solve(lambda, 1.4 + 0.02i, n_m, a, b, ...
        theta_inc, phi_inc, pol, opts);
    assert(lossy.info.validated);
    [Csca, Cext] = far_cross_sections(lossy);
    Cabs_surface = surface_absorption(lossy);
    [Cabs_volume_base, Cabs_volume] = volume_absorption(lossy);
    volume_quadrature_gap = relative_gap(Cabs_volume_base, Cabs_volume);
    assert(Csca >= 0 && Cext > Csca && Cabs_surface > 0 && ...
        Cabs_volume > 0);
    assert(volume_quadrature_gap <= quadrature_tolerance);
    assert(relative_gap(Cext, Csca + Cabs_surface) <= energy_tolerance);
    assert(relative_gap(Cext, Csca + Cabs_volume) <= energy_tolerance);
    assert(relative_gap(Cabs_surface, Cabs_volume) <= energy_tolerance);
    fprintf(['T7 %s lossy: Cext %.12g, Csca %.12g, ', ...
        'Cabs(P) %.12g, Cabs(V base) %.12g, Cabs(V refined) %.12g, ', ...
        'quadrature rel %.3g\n'], type, Cext, Csca, Cabs_surface, ...
        Cabs_volume_base, Cabs_volume, volume_quadrature_gap);
    solutions{k} = lossy;
end
end

function calibrate_far_field_with_mie
x = 5;
a = 1;
n_m = 1;
n_p = 1.5;
lambda = 2*pi*n_m*a/x;
theta_inc = 0.6;
phi_inc = -0.4;
pol = [1, 0.25i];
k_m = 2*pi*n_m/lambda;
[~, ~, khat, E0] = incident_plane_wave(k_m, n_m, theta_inc, ...
    phi_inc, pol, 0, 0, 0);
radius = 1e7*a;

[rhat, weights, grid_size] = spherical_grid(80, 160);
points = radius*rhat;
[Escattered, ~] = mie_field(lambda, n_p, n_m, a, theta_inc, ...
    phi_inc, pol, points(:, 1), points(:, 2), points(:, 3), ...
    'scattered');
intensity = reshape(row_norm(Escattered).^2, grid_size);
Csca = radius^2*weights.dphi*sum(weights.eta.*sum(intensity, 2))/ ...
    norm(E0)^2;

forward = radius*khat;
[Eforward, ~] = mie_field(lambda, n_p, n_m, a, theta_inc, ...
    phi_inc, pol, forward(1), forward(2), forward(3), 'scattered');
amplitude = radius*exp(-1i*k_m*radius)*Eforward;
Cext = (4*pi/k_m)*imag(sum(amplitude.*conj(E0)))/norm(E0)^2;

[Qext, Qsca] = mie_efficiencies(x, n_p/n_m);
assert_relative(Csca, pi*a^2*Qsca, 1e-6);
assert_relative(Cext, pi*a^2*Qext, 1e-6);
end

function [Csca, Cext] = far_cross_sections(sol)
eta_count = 2*sol.L + 12;
phi_count = 2*(2*sol.M + 5);
[rhat, weights, grid_size] = spherical_grid(eta_count, phi_count);
radius = 1e7*sol.semifocal;
points = radius*rhat;
[Escattered, ~] = spheroid_eval(sol, points(:, 1), points(:, 2), ...
    points(:, 3), 'scattered', 'auto');
intensity = reshape(row_norm(Escattered).^2, grid_size);
Csca = radius^2*weights.dphi*sum(weights.eta.*sum(intensity, 2))/ ...
    sol.E0_norm^2;

[~, ~, khat, E0] = incident_plane_wave(sol.km, sol.n_m, ...
    sol.theta_inc, sol.phi_inc, sol.pol, 0, 0, 0);
forward = radius*khat;
[Eforward, ~] = spheroid_eval(sol, forward(1), forward(2), ...
    forward(3), 'scattered', 'auto');
amplitude = radius*exp(-1i*sol.km*radius)*Eforward;
Cext = (4*pi/sol.km)*imag(sum(amplitude.*conj(E0)))/ ...
    sol.E0_norm^2;
end

function Cabs = surface_absorption(sol)
eta_count = 2*sol.L + 12;
phi_count = 2*(2*sol.M + 5);
[rhat, weights, grid_size] = spherical_grid(eta_count, phi_count);
radius = 1.25*max(sol.a, sol.b);
points = radius*rhat;
[E, H] = spheroid_eval(sol, points(:, 1), points(:, 2), ...
    points(:, 3), 'total', 'auto');
outward_flux = real(sum(cross(E, conj(H), 2).*rhat, 2));
outward_flux = reshape(outward_flux, grid_size);
integral = radius^2*weights.dphi* ...
    sum(weights.eta.*sum(outward_flux, 2));
Cabs = -integral/(sol.n_m*sol.E0_norm^2);
end

function [base, refined] = volume_absorption(sol)
if strcmp(sol.type, 'prolate')
    xi_minimum = 1;
else
    xi_minimum = 0;
end
xi_count = max(24, ceil(abs(sol.c_int)*(sol.xi0 - xi_minimum)) + 20);
eta_count = sol.L + 12;
phi_count = 2*(2*sol.M + 5);
base = volume_absorption_grid(sol, xi_minimum, xi_count, eta_count, ...
    phi_count);
% The phi grid exactly resolves Cartesian products through |Delta m| <= 2M+2.
refinement = 1.25;
refined = volume_absorption_grid(sol, xi_minimum, ...
    ceil(refinement*xi_count), ceil(refinement*eta_count), phi_count);
end

function Cabs = volume_absorption_grid(sol, xi_minimum, xi_count, ...
        eta_count, phi_count)
[xi_node, xi_weight] = gauss_legendre(xi_count);
[eta, eta_weight] = gauss_legendre(eta_count);
xi = xi_minimum + (sol.xi0 - xi_minimum)*(xi_node + 1)/2;
xi_weight = (sol.xi0 - xi_minimum)*xi_weight/2;
phi = 2*pi*(0:phi_count - 1)/phi_count;
[Xi, Eta, Phi] = ndgrid(xi, eta, phi);
if strcmp(sol.type, 'prolate')
    alpha = Xi.^2 - 1;
    D = Xi.^2 - Eta.^2;
else
    alpha = Xi.^2 + 1;
    D = Xi.^2 + Eta.^2;
end
rho = sqrt(max(0, alpha.*(1 - Eta.^2)));
X = sol.semifocal*rho.*cos(Phi);
Y = sol.semifocal*rho.*sin(Phi);
Z = sol.semifocal*Xi.*Eta;
[E, ~] = spheroid_eval(sol, X, Y, Z, 'total', 'auto');
intensity = reshape(row_norm(E).^2, size(Xi));
[Wxi, Weta, ~] = ndgrid(xi_weight, eta_weight, phi);
integrand = Wxi.*Weta.*D.*intensity;
field_integral = sol.semifocal^3*(2*pi/phi_count)*sum(integrand(:));
Cabs = sol.k0*imag(sol.n_p^2)*field_integral/ ...
    (sol.n_m*sol.E0_norm^2);
end

function test_boundary_traces(solutions)
for k = 1:numel(solutions)
    sol = solutions{k};
    eta_count = 2*(sol.L + 10);
    phi_count = 2*(2*sol.M + 5);
    eta = cos((2*(1:eta_count) - 1)*pi/(2*eta_count)).';
    phi = 2*pi*((0:phi_count - 1) + 0.5)/phi_count;
    [Eta, Phi] = ndgrid(eta, phi);
    rho = sol.b*sqrt(max(0, 1 - Eta.^2));
    X = rho.*cos(Phi);
    Y = rho.*sin(Phi);
    Z = sol.a*Eta;
    [Ein, Hin] = spheroid_eval(sol, X, Y, Z, 'total', 'in');
    [Eout, Hout] = spheroid_eval(sol, X, Y, Z, 'total', 'out');
    unit_normal = [X(:)/sol.b^2, Y(:)/sol.b^2, Z(:)/sol.a^2];
    unit_normal = unit_normal./row_norm(unit_normal);
    dE = Eout - Ein;
    dH = Hout - Hin;
    dE_tangent = dE - sum(dE.*unit_normal, 2).*unit_normal;
    dH_tangent = dH - sum(dH.*unit_normal, 2).*unit_normal;
    Eout_normal = sum(Eout.*unit_normal, 2);
    Ein_normal = sum(Ein.*unit_normal, 2);
    electric = max(row_norm(dE_tangent))/sol.E0_norm;
    magnetic = max(row_norm(dH_tangent))/ ...
        (sol.n_m*sol.E0_norm);
    normal_residual = max(abs(sol.n_m^2*Eout_normal - ...
        sol.n_p^2*Ein_normal))/ ...
        (sol.n_m^2*sol.E0_norm);
    assert(electric <= sol.opts.tol);
    assert(magnetic <= sol.opts.tol);
    assert(normal_residual <= sol.opts.tol);
    fprintf(['T8 %s: tangential E %.3g, tangential H %.3g, ', ...
        'normal n^2E %.3g\n'], sol.type, electric, magnetic, ...
        normal_residual);
end
end

function [rhat, weights, grid_size] = spherical_grid(eta_count, phi_count)
[eta, eta_weight] = gauss_legendre(eta_count);
phi = 2*pi*(0:phi_count - 1)/phi_count;
[Eta, Phi] = ndgrid(eta, phi);
sin_theta = sqrt(max(0, 1 - Eta.^2));
rhat = [sin_theta(:).*cos(Phi(:)), sin_theta(:).*sin(Phi(:)), Eta(:)];
weights = struct('eta', eta_weight, 'dphi', 2*pi/phi_count);
grid_size = size(Eta);
end

function [a, b] = axes_for_type(type, maximum_axis)
if strcmp(type, 'prolate')
    a = maximum_axis;
    b = maximum_axis/2;
else
    a = maximum_axis/2;
    b = maximum_axis;
end
end

function [Lx, Lz] = depolarization_factors(type, a, b)
if strcmp(type, 'prolate')
    eccentricity_squared = 1 - (b/a)^2;
    eccentricity = sqrt(eccentricity_squared);
    Lz = ((1 - eccentricity_squared)/eccentricity_squared)* ...
        (log((1 + eccentricity)/(1 - eccentricity))/(2*eccentricity) - 1);
    Lx = (1 - Lz)/2;
else
    eccentricity_squared = 1 - (a/b)^2;
    g = sqrt((1 - eccentricity_squared)/eccentricity_squared);
    Lx = (g/(2*eccentricity_squared))*(pi/2 - atan(g)) - g^2/2;
    Lz = 1 - 2*Lx;
end
end

function values = row_norm(array)
values = sqrt(sum(abs(array).^2, 2));
end

function value = relative_gap(first, second)
value = abs(first - second)/max([abs(first), abs(second), realmin]);
end

function assert_relative(actual, expected, tolerance)
scale = max([norm(actual(:)), norm(expected(:)), realmin]);
assert(norm(actual(:) - expected(:))/scale <= tolerance);
end

function assert_max_field(actual, expected, scale, tolerance)
assert(max(row_norm(actual - expected))/max(scale, realmin) <= tolerance);
end
