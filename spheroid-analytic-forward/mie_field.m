function [E, H] = mie_field(lambda, n_p, n_m, a, theta_inc, phi_inc, ...
    pol, X, Y, Z, field, side)
%MIE_FIELD Electromagnetic field of a homogeneous sphere.

if nargin < 11 || isempty(field)
    field = 'total';
end
if nargin < 12 || isempty(side)
    side = 'auto';
end
validateattributes(lambda, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'lambda');
validateattributes(n_p, {'numeric'}, {'finite', 'scalar'}, mfilename, 'n_p');
validateattributes(n_m, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'n_m');
validateattributes(a, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'a');
validateattributes(theta_inc, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'theta_inc');
validateattributes(phi_inc, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'phi_inc');
validateattributes(pol, {'numeric'}, ...
    {'finite', 'vector', 'numel', 2}, mfilename, 'pol');
validateattributes(X, {'numeric'}, {'real', 'finite'}, mfilename, 'X');
validateattributes(Y, {'numeric'}, {'real', 'finite'}, mfilename, 'Y');
validateattributes(Z, {'numeric'}, {'real', 'finite'}, mfilename, 'Z');
if imag(n_p) < 0 || n_p == 0
    error('mie_field:InvalidRefractiveIndex', ...
        'n_p must be nonzero with nonnegative imaginary part.');
end
if ~isequal(size(X), size(Y), size(Z))
    error('mie_field:InputSizeMismatch', 'X, Y, and Z must have equal sizes.');
end
field = option(field, {'total', 'scattered', 'incident'}, 'field');
side = option(side, {'auto', 'in', 'out'}, 'side');

k_m = 2*pi*n_m/lambda;
[E_inc, H_inc, khat] = incident_plane_wave( ...
    k_m, n_m, theta_inc, phi_inc, pol, X, Y, Z);
points = [X(:), Y(:), Z(:)];
r = hypot(hypot(points(:, 1), points(:, 2)), points(:, 3));
if strcmp(side, 'auto')
    outside = r/a >= 1 - 64*eps;
else
    if any(abs(r - a) > 1e-8*a)
        error('mie_field:InvalidForcedSide', ...
            'Forced side is allowed only within 1e-8*a of the boundary.');
    end
    outside = repmat(strcmp(side, 'out'), size(r));
end
if strcmp(field, 'incident')
    E = E_inc;
    H = H_inc;
    return
end

sin_theta = sin(theta_inc);
cos_theta = cos(theta_inc);
sin_phi = sin(phi_inc);
cos_phi = cos(phi_inc);
e_te = [-sin_phi, cos_phi, 0];
e_tm = [cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta];
e1 = cos_phi*e_tm - sin_phi*e_te;
e2 = sin_phi*e_tm + cos_phi*e_te;
R = [e1(:), e2(:), khat(:)];
points_local = points*R;
pol = pol(:).';

[~, ~, ~, coeff] = mie_efficiencies(k_m*a, n_p/n_m);
E_solution = complex(zeros(numel(r), 3));
H_solution = E_solution;
if any(outside)
    [E_local, H_local] = sphere_series(points_local(outside, :), k_m, ...
        n_m, pol, coeff, 'out');
    E_solution(outside, :) = E_local*R.';
    H_solution(outside, :) = H_local*R.';
end
inside = ~outside;
if any(inside)
    at_origin = inside & r == 0;
    regular = inside & ~at_origin;
    if any(regular)
        [E_local, H_local] = sphere_series(points_local(regular, :), ...
            2*pi*n_p/lambda, n_p, pol, coeff, 'in');
        E_solution(regular, :) = E_local*R.';
        H_solution(regular, :) = H_local*R.';
    end
    if any(at_origin)
        E_local = coeff.d(1)*[pol, 0];
        H_local = n_p*coeff.c(1)*[-pol(2), pol(1), 0];
        E_solution(at_origin, :) = repmat(E_local*R.', sum(at_origin), 1);
        H_solution(at_origin, :) = repmat(H_local*R.', sum(at_origin), 1);
    end
end

if strcmp(field, 'total')
    E = E_solution;
    H = H_solution;
    E(outside, :) = E(outside, :) + E_inc(outside, :);
    H(outside, :) = H(outside, :) + H_inc(outside, :);
else
    E = E_solution;
    H = H_solution;
    E(inside, :) = E(inside, :) - E_inc(inside, :);
    H(inside, :) = H(inside, :) - H_inc(inside, :);
end
end

function value = option(value, choices, name)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || ~ismember(value, choices)
    error('mie_field:InvalidOption', '%s must be %s.', name, ...
        strjoin(choices, ', '));
end
end

function [E, H] = sphere_series(points, k, n_medium, pol, coeff, region)
r = hypot(hypot(points(:, 1), points(:, 2)), points(:, 3));
rho = k*r;
mu = points(:, 3)./r;
sin_theta = hypot(points(:, 1), points(:, 2))./r;
phi = atan2(points(:, 2), points(:, 1));
phi(sin_theta == 0) = 0;
cos_phi = cos(phi);
sin_phi = sin(phi);
er = [sin_theta.*cos_phi, sin_theta.*sin_phi, mu];
etheta = [mu.*cos_phi, mu.*sin_phi, -sin_theta];
ephi = [-sin_phi, cos_phi, zeros(size(phi))];

p = pol;
q = [-pol(2), pol(1)];
Ap = p(1)*cos_phi + p(2)*sin_phi;
Bp = -p(1)*sin_phi + p(2)*cos_phi;
Aq = q(1)*cos_phi + q(2)*sin_phi;
Bq = -q(1)*sin_phi + q(2)*cos_phi;
Er = complex(zeros(size(r)));
Et = Er;
Ep = Er;
Hr = Er;
Ht = Er;
Hp = Er;
pi_nm1 = zeros(size(r));
pi_n = ones(size(r));

for j = 1:numel(coeff.n)
    n = coeff.n(j);
    tau = n*mu.*pi_n - (n + 1)*pi_nm1;
    [z, t, radial] = radial_terms(n, rho, region);
    factor = 1i^n*(2*n + 1)/(n*(n + 1));
    if strcmp(region, 'out')
        em = -coeff.b(j);
        en = 1i*coeff.a(j);
        hm = -n_medium*coeff.a(j);
        hn = 1i*n_medium*coeff.b(j);
    else
        em = coeff.c(j);
        en = -1i*coeff.d(j);
        hm = n_medium*coeff.d(j);
        hn = -1i*n_medium*coeff.c(j);
    end
    Er = Er + factor*en*(Ap.*sin_theta.*pi_n.*radial);
    Et = Et + factor*(em*Ap.*pi_n.*z + en*Ap.*tau.*t);
    Ep = Ep + factor*(em*Bp.*tau.*z + en*Bp.*pi_n.*t);
    Hr = Hr + factor*hn*(Aq.*sin_theta.*pi_n.*radial);
    Ht = Ht + factor*(hm*Aq.*pi_n.*z + hn*Aq.*tau.*t);
    Hp = Hp + factor*(hm*Bq.*tau.*z + hn*Bq.*pi_n.*t);
    pi_np1 = ((2*n + 1)/n)*mu.*pi_n - ((n + 1)/n)*pi_nm1;
    pi_nm1 = pi_n;
    pi_n = pi_np1;
end
E = Er.*er + Et.*etheta + Ep.*ephi;
H = Hr.*er + Ht.*etheta + Hp.*ephi;
end

function [z, t, q] = radial_terms(n, rho, region)
scale = sqrt(pi*rho/2);
psi_n = scale.*besselj(n + 0.5, rho);
psi_nm1 = scale.*besselj(n - 0.5, rho);
if strcmp(region, 'out')
    psi_n = psi_n + 1i*scale.*bessely(n + 0.5, rho);
    psi_nm1 = psi_nm1 + 1i*scale.*bessely(n - 0.5, rho);
end
z = psi_n./rho;
t = (psi_nm1 - n*psi_n./rho)./rho;
q = n*(n + 1)*psi_n./rho.^2;
end
