function test_scalar_approximation
% Characterize the fixed-polarization scalar approximation in the far field.
% This tests polarization leakage, not scalar-model amplitude or phase error.

lambda = 0.532; % um, vacuum wavelength
n_m = 1.335381534; % water at 20 deg C
n_p = 1.37; % homogeneous live-cell baseline; absorption neglected
opts = struct('tol', 1e-6, 'on_fail', 'error');
shapes = {'prolate', 'oblate'};
axes_ab = [5, 2.5; 2.5, 5]; % um
polarizations = {'x', 'y', 'z'};
theta_inc = [0, 0, pi/2];
pol = [1, 0; 0, 1; -1, 0];
component = [1, 2, 3];
narrow_offsets = linspace(-0.08, 0.08, 8);
wide_offsets = linspace(-1, 1, 10);
scalar_limit = 0.1;

for shape_index = 1:numel(shapes)
    a = axes_ab(shape_index, 1);
    b = axes_ab(shape_index, 2);
    radius = 1e4*max(a, b);
    for polarization_index = 1:numel(polarizations)
        sol = spheroid_solve(lambda, n_p, n_m, a, b, ...
            theta_inc(polarization_index), 0, ...
            pol(polarization_index, :), opts);
        assert(sol.info.validated, ...
            '%s %s-polarized solution was not validated.', ...
            shapes{shape_index}, polarizations{polarization_index});

        [~, ~, khat, E0] = incident_plane_wave(sol.km, n_m, ...
            theta_inc(polarization_index), 0, ...
            pol(polarization_index, :), 0, 0, 0);
        e0 = E0/norm(E0);
        expected = zeros(1, 3);
        expected(component(polarization_index)) = 1;
        assert(norm(e0 - expected) < 16*eps, ...
            'The requested polarization is not the intended Cartesian axis.');

        [narrow_leakage, narrow_retained] = polarization_fractions( ...
            sol, khat, e0, narrow_offsets, radius);
        [wide_leakage, wide_retained] = polarization_fractions( ...
            sol, khat, e0, wide_offsets, radius);
        fprintf(['scalar %-7s %s-pol: narrow leakage %.6f, retained %.6f; ', ...
            'wide leakage %.6f, retained %.6f\n'], ...
            shapes{shape_index}, polarizations{polarization_index}, ...
            narrow_leakage, narrow_retained, wide_leakage, wide_retained);
        assert(narrow_leakage < scalar_limit, ...
            ['%s %s polarization has %.3g discarded-vector leakage in ', ...
             'the narrow angular patch; scalar limit is %.3g.'], ...
            shapes{shape_index}, polarizations{polarization_index}, ...
            narrow_leakage, scalar_limit);
        assert(wide_leakage > scalar_limit, ...
            ['%s %s wide-patch leakage %.3g did not expose the expected ', ...
             'non-scalar vector field.'], ...
            shapes{shape_index}, polarizations{polarization_index}, ...
            wide_leakage);
    end
end

fprintf(['test_scalar_approximation: PASS ', ...
    '(fixed-polarization reduction only on the tested narrow patch)\n']);
end

function [leakage, retained] = polarization_fractions( ...
    sol, khat, e0, offsets, radius)
[first_offset, second_offset] = ndgrid(offsets, offsets);
second_transverse = cross(khat, e0);
directions = khat + first_offset(:)*e0 + ...
    second_offset(:)*second_transverse;
directions = directions./sqrt(sum(abs(directions).^2, 2));
points = radius*directions;
E = spheroid_eval(sol, points(:, 1), points(:, 2), points(:, 3), ...
    'scattered', 'auto');
field_norm = norm(E, 'fro');
assert(field_norm > realmin, 'The sampled scattered field is numerically zero.');
coefficient = E*conj(e0(:));
co_polarized = coefficient*e0;
leakage = norm(E - co_polarized, 'fro')/field_norm;
retained = norm(co_polarized, 'fro')/field_norm;
end
