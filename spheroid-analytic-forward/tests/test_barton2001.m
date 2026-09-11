function test_barton2001
% Far-field maxima from Barton, Appl. Opt. 40, 3598-3607 (2001).
% Source: ref/mie-scatt/ao-40-21-3598.pdf. Geometry and illumination are
% on p. 3602, Eq. (50) is on p. 3603, and the reference maxima are in
% the captions of Figs. 7 and 13 on pp. 3604-3605.

f = 1;
n_m = 1;
theta_inc = pi/6;
phi_inc = 0;
pol = [0, 1];
step_deg = 0.25;
psi_deg = (-180:step_deg:(180 - step_deg)).';
psi = (pi/180)*psi_deg;
radius = 1e7*f;
X = radius*sin(psi);
Y = zeros(size(psi));
Z = radius*cos(psi);

cases = { ...
    'prolate, n=1.33', 1.154701, sqrt(1.154701^2 - 1), ...
        27.494593, 1.33, 12.18744; ...
    'prolate, n=1.50', 1.154701, sqrt(1.154701^2 - 1), ...
        27.494593, 1.50, 14.48922; ...
    'oblate, n=1.33', 0.577350, sqrt(0.577350^2 + 1), ...
        21.822473, 1.33, 60.93952; ...
    'oblate, n=1.50', 0.577350, sqrt(0.577350^2 + 1), ...
        21.822473, 1.50, 51.01663};

for k = 1:size(cases, 1)
    name = cases{k, 1};
    a = cases{k, 2};
    b = cases{k, 3};
    h_ext = cases{k, 4};
    n_p = cases{k, 5};
    expected_peak = cases{k, 6};
    lambda = 2*pi/h_ext;

    sol = spheroid_solve(lambda, n_p, n_m, a, b, theta_inc, ...
        phi_inc, pol);
    assert(sol.info.validated, '%s did not produce a validated solution.', ...
        name);

    E_s = spheroid_eval(sol, X, Y, Z, 'scattered');
    S_r = radius^2*n_m*sum(abs(E_s).^2, 2) / ...
        (pi*f^2*sol.E0_norm^2);
    [peak, peak_index] = max(S_r);
    peak_angle_deg = psi_deg(peak_index);
    relative_error = abs(peak - expected_peak)/expected_peak;
    assert(relative_error < 1e-4, ...
        '%s: S_r,max %.8g differs from %.8g by %.3g.', ...
        name, peak, expected_peak, relative_error);

    % The caption values are sampled-plot maxima; refining the continuous
    % maximum changes the n=1.50 prolate value beyond the quoted digits.
    fprintf(['Barton %-18s S_r,max=%11.7f reference=%11.7f ', ...
        'relerr=%.3g angle=%8.3f deg\n'], name, peak, expected_peak, ...
        relative_error, peak_angle_deg);
end

fprintf('test_barton2001: PASS\n');
end
