function test_incident_plane_wave
% Focused plane-wave convention, curl, array, and validation checks.

k_m = 2;
n_m = 1.5;
Z = [0; pi/4];
[E, H, khat, E0] = incident_plane_wave( ...
    k_m, n_m, 0, 1.2, [1, 0], zeros(2, 1), zeros(2, 1), Z);
assert_close(khat, [0, 0, 1], 5e-15);
assert_close(E0, [1, 0, 0], 5e-15);
assert_close(E, [1, 0, 0; 1i, 0, 0], 5e-15);
assert_close(H, [0, n_m, 0; 0, 1i*n_m, 0], 5e-15);

theta = 0.8;
phi = -0.3;
pol = [1 + 2i, -0.5i];
X = [0.2, -0.4, 0.1];
Y = [-0.1, 0.3, 0.8];
Z = [0.7, -0.2, 0.5];
[E, H, khat, E0] = incident_plane_wave( ...
    3.1, 1.2, theta, phi, pol, X, Y, Z);
expected_khat = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
e_te = [-sin(phi), cos(phi), 0];
e_tm = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)];
e1 = cos(phi)*e_tm - sin(phi)*e_te;
e2 = sin(phi)*e_tm + cos(phi)*e_te;
expected_E0 = pol(1)*e1 + pol(2)*e2;
expected_phase = exp(1i*3.1*(expected_khat(1)*X(:) + ...
    expected_khat(2)*Y(:) + expected_khat(3)*Z(:)));
expected_H0 = 1.2*cross(expected_khat, expected_E0);
assert_close(khat, expected_khat, 5e-15);
assert_close(E0, expected_E0, 5e-15);
assert_close(E, expected_phase .* expected_E0, 5e-14);
assert_close(H, expected_phase .* expected_H0, 5e-14);
assert(isequal(size(E), [3, 3]));
assert(isequal(size(H), [3, 3]));
assert(abs(khat*E0.') < 5e-15*max(1, norm(E0)));

check_curl_relation;

[E, H] = incident_plane_wave(1, 1, 0.4, 0.2, [0; 0], 1, 2, 3);
assert(isequal(E, zeros(1, 3)));
assert(isequal(H, zeros(1, 3)));

must_error(@() incident_plane_wave(0, 1, 0, 0, [1, 0], 0, 0, 0));
must_error(@() incident_plane_wave(1i, 1, 0, 0, [1, 0], 0, 0, 0));
must_error(@() incident_plane_wave(1, 0, 0, 0, [1, 0], 0, 0, 0));
must_error(@() incident_plane_wave(1, 1i, 0, 0, [1, 0], 0, 0, 0));
must_error(@() incident_plane_wave(1, 1, [0, 1], 0, [1, 0], 0, 0, 0));
must_error(@() incident_plane_wave(1, 1, 0, NaN, [1, 0], 0, 0, 0));
must_error(@() incident_plane_wave(1, 1, 0, 0, [1, 0, 0], 0, 0, 0));
must_error(@() incident_plane_wave(1, 1, 0, 0, ones(2), 0, 0, 0));
must_error(@() incident_plane_wave(1, 1, 0, 0, [1, NaN], 0, 0, 0));
must_error(@() incident_plane_wave(1, 1, 0, 0, [1, 0], [0, 1], [0; 1], [0, 1]));
must_error(@() incident_plane_wave(1, 1, 0, 0, [1, 0], 1i, 0, 0));
must_error(@() incident_plane_wave(1, 1, 0, 0, [1, 0], 0, Inf, 0));

fprintf('test_incident_plane_wave: PASS\n');
end

function check_curl_relation
k_m = 2.3;
n_m = 1.4;
theta = 0.6;
phi = 1.1;
pol = [0.7 - 0.2i, -0.4 + 0.3i];
p = [0.2, -0.3, 0.4];
h = 1e-5;

[Epx, ~] = incident_plane_wave(k_m, n_m, theta, phi, pol, p(1)+h, p(2), p(3));
[Emx, ~] = incident_plane_wave(k_m, n_m, theta, phi, pol, p(1)-h, p(2), p(3));
[Epy, ~] = incident_plane_wave(k_m, n_m, theta, phi, pol, p(1), p(2)+h, p(3));
[Emy, ~] = incident_plane_wave(k_m, n_m, theta, phi, pol, p(1), p(2)-h, p(3));
[Epz, ~] = incident_plane_wave(k_m, n_m, theta, phi, pol, p(1), p(2), p(3)+h);
[Emz, ~] = incident_plane_wave(k_m, n_m, theta, phi, pol, p(1), p(2), p(3)-h);
[~, H] = incident_plane_wave(k_m, n_m, theta, phi, pol, p(1), p(2), p(3));

dE_dx = (Epx - Emx)/(2*h);
dE_dy = (Epy - Emy)/(2*h);
dE_dz = (Epz - Emz)/(2*h);
curl_E = [dE_dy(3) - dE_dz(2), ...
          dE_dz(1) - dE_dx(3), ...
          dE_dx(2) - dE_dy(1)];
expected = 1i*(k_m/n_m)*H;
assert(norm(curl_E - expected)/norm(expected) < 1e-9);
end

function assert_close(actual, expected, tolerance)
scale = max(1, max(abs(expected(:))));
assert(max(abs(actual(:) - expected(:))) <= tolerance*scale);
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
