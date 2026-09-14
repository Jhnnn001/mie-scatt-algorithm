function test_solution_validation
% Reject inconsistent active-sector and outgoing radial data.

sol = minimal_solution();
[E, H] = spheroid_eval(sol, 0, 0, 0, 'incident');
assert(all(isfinite([E(:); H(:)])));

bad = sol;
bad.mu_data{1}.mu = 1;
must_error_id(@() spheroid_eval(bad, 0, 0, 0, 'incident'), ...
    'spheroid_eval:InvalidSolution');

bad = sol;
bad.mu_data{1}.rod_ext = rmfield(bad.mu_data{1}.rod_ext, 'solution');
must_error_id(@() spheroid_eval(bad, 0, 0, 0, 'incident'), ...
    'spheroid_eval:InvalidSolution');

bad = sol;
bad.mu_data{1}.rod_ext.solution = [];
must_error_id(@() spheroid_eval(bad, 0, 0, 0, 'incident'), ...
    'spheroid_eval:InvalidSolution');

bad = sol;
bad.mu_data{1}.rod_ext.qtilde_xi0 = [];
must_error_id(@() spheroid_eval(bad, 0, 0, 0, 'incident'), ...
    'spheroid_eval:InvalidSolution');

bad = sol;
bad.mu_data{1}.rod_ext.asymptotic_coefficients(1) = NaN;
must_error_id(@() spheroid_eval(bad, 0, 0, 0, 'incident'), ...
    'spheroid_eval:InvalidSolution');

bad = sol;
bad.mu_data{1}.rod_ext.arbitrary_scale = false;
must_error_id(@() spheroid_eval(bad, 0, 0, 0, 'incident'), ...
    'spheroid_eval:InvalidSolution');

bad = sol;
bad.mu_data{1}.rod_ext.normalization_xi = sol.xi0 + 1;
must_error_id(@() spheroid_eval(bad, 0, 0, 0, 'incident'), ...
    'spheroid_eval:InvalidSolution');
fprintf('test_solution_validation: PASS\n');
end

function sol = minimal_solution()
type = 'prolate';
mu = 0;
L = 0;
semifocal = 1;
xi0 = 2;
a = semifocal*xi0;
b = semifocal*sqrt(xi0^2 - 1);
c_ext = 0.2;
c_int = 0.24;
eg_ext = sph_eigen(type, mu, c_ext, L);
eg_int = sph_eigen(type, mu, c_int, L);
[rod_ext, check] = sph_radial_ode(eg_ext, c_ext, xi0, 1.5);
assert(check.ok);

mode = struct('m', 0, 'mu', 0, 'a', 0, 'b', 0, 'c', 0, 'd', 0, ...
    'analytically_zeroed', false);
data = struct('mu', mu, 'eg_ext', eg_ext, 'eg_int', eg_int, ...
    'rod_ext', rod_ext);
sol = struct('route', 'spheroid', 'type', type, 'a', a, 'b', b, ...
    'semifocal', semifocal, 'xi0', xi0, 'c_ext', c_ext, ...
    'c_int', c_int, 'km', c_ext/semifocal, 'n_m', 1, 'n_p', 1.2, ...
    'theta_inc', 0, 'phi_inc', 0, 'pol', [1, 0], 'L', L, 'M', 0, ...
    'modes', mode, 'mu_data', {{data}});
end

function must_error_id(f, identifier)
try
    f();
catch exception
    assert(strcmp(exception.identifier, identifier), ...
        'Expected %s, got %s.', identifier, exception.identifier);
    return
end
error('test_solution_validation:MissingError', ...
    'Expected error %s.', identifier);
end
