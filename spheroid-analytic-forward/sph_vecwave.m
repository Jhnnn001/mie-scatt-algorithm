function [M, N] = sph_vecwave(W, V, g, m, c)
%SPH_VECWAVE Phase-free vector spheroidal wavefunctions in coordinate bases.

narginchk(5, 5);
angular_names = {'W0', 'W1', 'W2', 'Wm'};
radial_names = {'V0', 'V1', 'V2', 'Va', 'Vm'};
if ~isstruct(W) || ~all(isfield(W, angular_names))
    error('sph_vecwave:InvalidAngularData', ...
        'W must be a structure returned by sph_angular.');
end
if ~isstruct(V) || ~all(isfield(V, radial_names))
    error('sph_vecwave:InvalidRadialData', ...
        'V must be a structure returned by sph_radial.');
end
coordinate_names = {'eta', 'xi', 'phi', 'sigma', 'alpha', 'beta', 'D', ...
    'A_eta_prime', 'B_prime', 'D_eta', 'D_xi'};
if ~isstruct(g) || ~all(isfield(g, coordinate_names))
    error('sph_vecwave:InvalidCoordinates', ...
        'g must be a structure returned by sph_coords.');
end
if ~(isnumeric(m) && isreal(m) && isfinite(m) && isscalar(m) && ...
        m == fix(m) && abs(m) <= 300)
    error('sph_vecwave:InvalidAzimuthalIndex', ...
        'm must be an integer scalar with abs(m) <= 300.');
end
if ~(isnumeric(c) && isscalar(c) && isfinite(c) && c ~= 0)
    error('sph_vecwave:InvalidParameter', ...
        'c must be a finite nonzero numeric scalar.');
end

shape = size(W.W0);
if numel(shape) > 2 || any(shape == 0)
    error('sph_vecwave:SizeMismatch', ...
        'Angular data must be a nonempty modes-by-points matrix.');
end
for name = angular_names
    validate_modal_array(W.(name{1}), shape);
end
for name = radial_names
    validate_modal_array(V.(name{1}), shape);
end
points = shape(2);
for name = coordinate_names
    value = g.(name{1});
    if ~(isnumeric(value) && isvector(value) && numel(value) == points && ...
            isreal(value) && all(isfinite(value(:))))
        error('sph_vecwave:InvalidCoordinates', ...
            'Coordinate fields must be finite real vectors matching the point count.');
    end
end

eta = g.eta(:).';
xi = g.xi(:).';
alpha = g.alpha(:).';
beta = g.beta(:).';
D = g.D(:).';
A_prime = g.A_eta_prime(:).';
B_prime = g.B_prime(:).';
D_eta = g.D_eta(:).';
D_xi = g.D_xi(:).';
sigma_values = g.sigma(:).';
if all(sigma_values == 1)
    sigma = 1;
elseif all(sigma_values == -1)
    sigma = -1;
else
    error('sph_vecwave:InvalidCoordinates', ...
        'Coordinate data do not identify one spheroidal type.');
end
coordinate_scale = max(1, max(abs([xi, alpha, D, B_prime])));
coordinate_tolerance = 1e-11*coordinate_scale;
if any(abs(eta) > 1) || any(xi < max(0, sigma)) || ...
        any(alpha < (1 - sigma)/2) || ...
        any(beta < 0 | beta > 1) || any(D <= 0) || ...
        max(abs(beta - (1 - eta.^2))) > coordinate_tolerance || ...
        max(abs(alpha - (xi.^2 - sigma))) > coordinate_tolerance || ...
        max(abs(D - (alpha + sigma*beta))) > coordinate_tolerance || ...
        max(abs(A_prime - sigma*(1 - 3*eta.^2))) > coordinate_tolerance || ...
        max(abs(B_prime - (3*xi.^2 - sigma))) > coordinate_tolerance || ...
        max(abs(D_eta + 2*sigma*eta)) > coordinate_tolerance || ...
        max(abs(D_xi - 2*xi)) > coordinate_tolerance
    error('sph_vecwave:InvalidCoordinates', ...
        'Coordinate values and metric fields are inconsistent.');
end

sqrt_alpha = sqrt(alpha);
sqrt_beta = sqrt(beta);
sqrt_D = sqrt(D);
P00 = W.W0.*V.V0;
P10 = W.W1.*V.V0;
P01 = W.W0.*V.V1;
P20 = W.W2.*V.V0;
P02 = W.W0.*V.V2;
P11 = W.W1.*V.V1;
PmA = W.Wm.*V.Va;

s = sign(m);
M = struct( ...
    'eta', -1i*s*PmA.*(xi.*sqrt(alpha./D)), ...
    'xi', 1i*s*PmA.*(sigma*eta.*sqrt(beta./D)), ...
    'phi', (xi.*sqrt_alpha.*P10 - ...
        sigma*eta.*sqrt_beta.*P01)./D);

eta_euler = sigma*eta.*beta.*P10 + ...
    xi.*sqrt(alpha.*beta).*P01;
xi_euler = sigma*eta.*sqrt(alpha.*beta).*P10 + ...
    xi.*alpha.*P01;
Neta = P10 + (A_prime.*P10 + sigma*eta.*P20 + ...
    xi.*sqrt_alpha.*P11)./D - D_eta.*eta_euler./D.^2;
Nxi = P01 + (B_prime.*P01 + xi.*P02 + ...
    sigma*eta.*sqrt_beta.*P11)./D - D_xi.*xi_euler./D.^2;
N = struct( ...
    'eta', Neta./(c*sqrt_D) + ...
        c*sigma*eta.*sqrt(beta./D).*P00, ...
    'xi', Nxi./(c*sqrt_D) + c*xi.*sqrt(alpha./D).*P00, ...
    'phi', (1i*s/c)*(PmA + ...
        (sigma*eta.*W.W1.*V.Vm + xi.*W.Wm.*V.V1)./D));
end

function validate_modal_array(value, shape)
if ~(isnumeric(value) && isequal(size(value), shape) && ...
        all(isfinite(real(value(:)))) && all(isfinite(imag(value(:)))))
    error('sph_vecwave:SizeMismatch', ...
        'All angular and radial fields must have one common finite size.');
end
end
