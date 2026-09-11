function test_sph_coords
% Focused coordinate, basis, singular-set, and validation checks.

delta = 1e-8;

for scale = [1e-200, 1, 1e200]
    actual = sph_semifocal(scale, scale/2);
    expected = scale*sqrt(3)/2;
    assert(abs(actual/expected - 1) < 4*eps);
end
assert(sph_semifocal(1, 1) == 0);
assert(abs(sph_semifocal(uint8(1), uint8(2)) - sqrt(3)) < 4*eps);
must_error(@() sph_semifocal(0, 1));

g = sph_coords('prolate', [0, 0, 0], [0, 0, 0], [-0.5, 0, 0.5]);
assert_close(g.xi, ones(3, 1), 0);
assert_close(g.eta, [-0.5; 0; 0.5], 0);
assert_close(g.alpha, zeros(3, 1), 0);
assert_close(g.beta, [0.75; 1; 0.75], 0);
assert_close(g.D, g.beta, 0);
assert(~any(g.nudged));
assert(all(isinf(g.h_xi)));
check_basis(g, 3);

g = sph_coords('prolate', 1e-200, 0, 0.5);
assert_close(g.xi, 1, 0);
assert_close(g.eta, 0.5, 0);
assert_close(g.alpha, 0, 0);
assert_close(g.beta, 0.75, 0);
assert_close(g.D, 0.75, 0);
check_basis(g, 1);

g = sph_coords('prolate', [0; 0], [0; 0], [-2; 2]);
assert(~any(g.nudged));
assert_close(g.xi, [2; 2], 0);
assert_close(g.eta, [-1; 1], 0);
assert_close(g.beta, [0; 0], 0);
check_basis(g, 2);

g = sph_coords('prolate', [0; 0], [0; 0], [-1; 1]);
assert(all(g.nudged));
assert(all(g.D > 0));
check_basis(g, 2);

g = sph_coords('prolate', delta/2, 0, 1);
assert(g.nudged);
assert(g.D > 0);
check_basis(g, 1);

g = sph_coords('prolate', delta, 0, 1);
assert(~g.nudged);

near_delta = delta - eps(delta);
x = [-near_delta; near_delta; -near_delta; near_delta];
g = sph_coords('prolate', x, zeros(4, 1), [1; 1; -1; -1]);
assert(all(g.nudged));
assert(max(abs(sqrt(g.alpha .* g.beta)/delta - 1)) < 1e-12);
assert(all(g.D >= delta));
assert_close(cos(g.phi), sign(x), 5e-15);
check_basis(g, 4);

phi = 0.7;
rho = 1.5;
g = sph_coords('prolate', rho*cos(phi), rho*sin(phi), 1);
assert_close(g.xi, 2, 5e-15);
assert_close(g.eta, 0.5, 5e-15);
assert_close(g.phi, phi, 5e-15);
assert_close(g.alpha, 3, 5e-15);
assert_close(g.beta, 0.75, 5e-15);
assert_close(g.D, 3.75, 5e-15);
assert_close(g.A_eta, 0.375, 5e-15);
assert_close(g.A_eta_prime, 0.25, 5e-15);
assert_close(g.B, 6, 5e-15);
assert_close(g.B_prime, 11, 5e-15);
assert_close(g.D_eta, -1, 5e-15);
assert_close(g.D_xi, 4, 5e-15);
check_basis(g, 1);

g = sph_coords("oblate", [0; 0.6; -0.8], zeros(3, 1), zeros(3, 1));
assert_close(g.xi, zeros(3, 1), 0);
assert_close(g.eta, [1; 0.8; 0.6], 5e-15);
assert(~any(g.nudged));
check_basis(g, 3);

g = sph_coords('oblate', [0; 0; 0], [0; 0; 0], [-2; 0; 2]);
assert(~any(g.nudged));
assert_close(g.xi, [2; 0; 2], 5e-15);
assert_close(g.eta, [-1; 1; 1], 5e-15);
check_basis(g, 3);

tiny_z = [1e-200; -1e-200; 1e-200; -1e-200];
axis_z = linspace(-5, 9, 176).'/4;
g = sph_coords('oblate', zeros(size(axis_z)), zeros(size(axis_z)), axis_z);
assert(all(abs(g.eta) <= 1));
check_basis(g, numel(axis_z));

g = sph_coords('oblate', [0.5; 0.5; 2; 2], zeros(4, 1), tiny_z);
assert(all(g.xi(1:2) > 0));
assert(all(abs(g.eta(3:4)) > 0));
assert(max(abs(g.xi .* g.eta ./ tiny_z - 1)) < 5e-15);
check_basis(g, 4);

g = sph_coords('oblate', [1; 0; -1], [0; 1; 0], zeros(3, 1));
assert(all(g.nudged));
assert(all(g.D > 0));
check_basis(g, 3);

ring_phi = [0; pi/2; 0.73];
g = sph_coords('oblate', cos(ring_phi), sin(ring_phi), zeros(3, 1));
assert(all(g.nudged));
assert(all(isfinite(g.D)) && all(g.D > 0));
azimuth_error = atan2(sin(g.phi - ring_phi), cos(g.phi - ring_phi));
assert(max(abs(azimuth_error)) < 20*eps);
assert((max(g.D) - min(g.D))/mean(g.D) < 1e-7);
check_basis(g, 3);

ring_x = -delta/2;
ring_y = sqrt(1 - ring_x^2);
g = sph_coords('oblate', ring_x, ring_y, 0);
assert(g.nudged);
assert(g.D > 0);
check_basis(g, 1);

g = sph_coords('oblate', -delta, 1, 0);
assert(g.nudged);
assert(g.D > 0);
check_basis(g, 1);

g = sph_coords('oblate', 1, 0, delta/2);
assert(g.nudged);
assert(g.D > 0);
check_basis(g, 1);

g = sph_coords('oblate', 1, 0, delta);
assert(~g.nudged);

phi_nudge = 0.4;
rho_nudge = [1 - near_delta; 1 + near_delta];
g = sph_coords('oblate', rho_nudge*cos(phi_nudge), ...
    rho_nudge*sin(phi_nudge), zeros(2, 1));
assert(all(g.nudged));
rho_after = sqrt(g.alpha .* g.beta);
assert(max(abs(abs(rho_after - 1)/delta - 1)) < 1e-8);
assert(all(g.D >= delta));
assert_close(g.phi, repmat(phi_nudge, 2, 1), 5e-15);
check_basis(g, 2);

rho = sqrt(15)/2;
g = sph_coords('oblate', rho*cos(phi), rho*sin(phi), 1);
assert_close(g.xi, 2, 5e-15);
assert_close(g.eta, 0.5, 5e-15);
assert_close(g.phi, phi, 5e-15);
assert_close(g.alpha, 5, 5e-15);
assert_close(g.beta, 0.75, 5e-15);
assert_close(g.D, 4.25, 5e-15);
assert_close(g.A_eta, -0.375, 5e-15);
assert_close(g.A_eta_prime, -0.25, 5e-15);
assert_close(g.B, 10, 5e-15);
assert_close(g.B_prime, 13, 5e-15);
assert_close(g.D_eta, 1, 5e-15);
assert_close(g.D_xi, 4, 5e-15);
check_basis(g, 1);

r = 1e7;
theta = 0.8;
phi = -0.4;
x = r*sin(theta)*cos(phi);
y = r*sin(theta)*sin(phi);
z = r*cos(theta);
gp = sph_coords('prolate', x, y, z);
go = sph_coords('oblate', x, y, z);
assert(abs(gp.alpha/(r^2 - 1) - 1) < 1e-12);
assert(abs(go.alpha/(r^2 + 1) - 1) < 1e-12);
assert(abs(gp.beta/sin(theta)^2 - 1) < 1e-12);
assert(abs(go.beta/sin(theta)^2 - 1) < 1e-12);

must_error(@() sph_coords('sphere', 0, 0, 0));
must_error(@() sph_coords(1, 0, 0, 0));
must_error(@() sph_coords('prolate', [0, 1], [0; 1], [0, 1]));
must_error(@() sph_coords('prolate', 1i, 0, 0));
must_error(@() sph_coords('oblate', 0, NaN, 0));

fprintf('test_sph_coords: PASS\n');
end

function check_basis(g, n)
scalar_fields = {'eta', 'xi', 'phi', 'alpha', 'beta', 'D', ...
    'h_eta', 'h_xi', 'h_phi', 'A_eta', 'A_eta_prime', ...
    'B', 'B_prime', 'D_eta', 'D_xi', 'nudged'};
for k = 1:numel(scalar_fields)
    assert(isequal(size(g.(scalar_fields{k})), [n, 1]));
end
assert(isequal(size(g.e_eta), [n, 3]));
assert(isequal(size(g.e_xi), [n, 3]));
assert(isequal(size(g.e_phi), [n, 3]));
assert(all(isfinite(g.e_eta(:))));
assert(all(isfinite(g.e_xi(:))));
assert(all(isfinite(g.e_phi(:))));
assert(max(abs(sum(g.e_eta.^2, 2) - 1)) < 2e-14);
assert(max(abs(sum(g.e_xi.^2, 2) - 1)) < 2e-14);
assert(max(abs(sum(g.e_phi.^2, 2) - 1)) < 2e-14);
assert(max(abs(sum(g.e_eta .* g.e_xi, 2))) < 2e-14);
assert(max(abs(cross(g.e_eta, g.e_xi, 2) - g.e_phi), [], 'all') < 2e-14);
end

function assert_close(actual, expected, tolerance)
scale = max(1, max(abs(expected(:))));
assert(max(abs(actual(:) - expected(:))) <= tolerance * scale);
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
