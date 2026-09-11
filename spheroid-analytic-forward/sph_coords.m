function g = sph_coords(type, x, y, z)
%SPH_COORDS Stable prolate or oblate spheroidal coordinates in focal units.

if isstring(type) && isscalar(type)
    type = char(type);
end
if ~(ischar(type) && isrow(type) && ...
        (strcmp(type, 'prolate') || strcmp(type, 'oblate')))
    error('sph_coords:InvalidType', ...
        'type must be ''prolate'' or ''oblate''.');
end
validateattributes(x, {'numeric'}, {'real', 'finite'}, mfilename, 'x');
validateattributes(y, {'numeric'}, {'real', 'finite'}, mfilename, 'y');
validateattributes(z, {'numeric'}, {'real', 'finite'}, mfilename, 'z');
if ~isequal(size(x), size(y), size(z))
    error('sph_coords:InputSizeMismatch', ...
        'x, y, and z must have equal sizes.');
end

x = x(:);
y = y(:);
z = z(:);
rho_original = hypot(x, y);
delta = 1e-8;

if strcmp(type, 'prolate')
    sigma = ones(size(x));
    abs_z = abs(z);
    nudged = hypot(rho_original, abs_z - 1) < delta;
    dz = abs_z - 1;
    target_rho = zeros(size(x));
    target_rho(nudged) = sqrt(max(0, ...
        (delta - abs(dz(nudged))) .* (delta + abs(dz(nudged)))));
    with_azimuth = nudged & rho_original > 0;
    radial_scale = target_rho(with_azimuth) ./ rho_original(with_azimuth);
    x(with_azimuth) = x(with_azimuth) .* radial_scale;
    y(with_azimuth) = y(with_azimuth) .* radial_scale;
    on_axis = nudged & rho_original == 0;
    x(on_axis) = target_rho(on_axis);
    rho = hypot(x, y);

    r_near = hypot(rho, abs_z - 1);
    r_far = hypot(rho, abs_z + 1);
    xi_minus_one = zeros(size(x));
    inside = abs_z <= 1;
    xi_minus_one(inside) = rho(inside).^2 / 2 .* (...
        1 ./ (r_far(inside) + 1 + abs_z(inside)) + ...
        1 ./ (r_near(inside) + 1 - abs_z(inside)));
    outside = ~inside;
    xi_minus_one(outside) = abs_z(outside) - 1 + ...
        rho(outside).^2 / 2 .* (...
        1 ./ (r_far(outside) + abs_z(outside) + 1) + ...
        1 ./ (r_near(outside) + abs_z(outside) - 1));

    xi = 1 + xi_minus_one;
    alpha = xi_minus_one .* (xi + 1);
    segment = alpha == 0 & abs_z < 1;
    beta = zeros(size(x));
    regular = ~segment;
    beta(regular) = rho(regular).^2 ./ alpha(regular);
    beta(segment) = 1 - z(segment).^2;
    beta = min(1, max(0, beta));
    D = alpha + beta;

    eta = zeros(size(x));
    direct = beta >= 0.5;
    eta(direct) = sign(z(direct)) .* ...
        (r_far(direct) - r_near(direct)) / 2;
    stable = ~direct;
    one_minus_abs_eta = beta(stable) ./ ...
        (1 + sqrt(1 - beta(stable)));
    eta(stable) = sign(z(stable)) .* (1 - one_minus_abs_eta);
    eta(segment) = z(segment);

    A_eta = eta .* beta;
    A_eta_prime = 1 - 3*eta.^2;
    B_prime = 3*xi.^2 - 1;
    D_eta = -2*eta;
else
    sigma = -ones(size(x));
    nudged = hypot(rho_original - 1, z) < delta;
    dr = rho_original(nudged) - 1;
    side = ones(size(dr));
    side(dr < 0) = -1;
    target_dr = sqrt(max(0, ...
        (delta - abs(z(nudged))) .* (delta + abs(z(nudged)))));
    target_rho = 1 + side .* target_dr;
    radial_scale = target_rho ./ rho_original(nudged);
    x(nudged) = x(nudged) .* radial_scale;
    y(nudged) = y(nudged) .* radial_scale;
    rho = hypot(x, y);

    A = (rho - 1) .* (rho + 1) + z.^2;
    dominant_x = abs(x) >= abs(y);
    A_stable = zeros(size(A));
    A_stable(dominant_x) = (x(dominant_x) - 1) .* ...
        (x(dominant_x) + 1) + y(dominant_x).^2 + z(dominant_x).^2;
    A_stable(~dominant_x) = x(~dominant_x).^2 + ...
        (y(~dominant_x) - 1) .* (y(~dominant_x) + 1) + z(~dominant_x).^2;
    A(nudged) = A_stable(nudged);
    D = hypot(A, 2*z);

    positive_A = A >= 0;
    xi = zeros(size(x));
    eta_abs = zeros(size(x));
    xi(positive_A) = sqrt(max(0, ...
        (D(positive_A) + A(positive_A)) / 2));
    eta_abs(positive_A) = abs(z(positive_A)) ./ xi(positive_A);
    negative_A = ~positive_A;
    eta_abs(negative_A) = sqrt(max(0, ...
        (D(negative_A) - A(negative_A)) / 2));
    xi(negative_A) = abs(z(negative_A)) ./ eta_abs(negative_A);
    eta_abs = min(1, eta_abs); % Roundoff on the symmetry axis can exceed one.

    eta_sign = ones(size(z));
    eta_sign(z < 0) = -1;
    eta = eta_sign .* eta_abs;
    alpha = xi.^2 + 1;
    beta = min(1, max(0, rho.^2 ./ alpha));

    A_eta = -eta .* beta;
    A_eta_prime = -(1 - 3*eta.^2);
    B_prime = 3*xi.^2 + 1;
    D_eta = 2*eta;
end

phi = atan2(y, x);
sqrt_alpha_over_D = sqrt(alpha ./ D);
sqrt_beta_over_D = sqrt(beta ./ D);
cos_phi = cos(phi);
sin_phi = sin(phi);

g.eta = eta;
g.xi = xi;
g.phi = phi;
g.sigma = sigma;
g.alpha = alpha;
g.beta = beta;
g.D = D;
g.h_eta = sqrt(D ./ beta);
g.h_xi = sqrt(D ./ alpha);
g.h_phi = sqrt(alpha .* beta);
g.e_eta = [-eta .* sqrt_alpha_over_D .* cos_phi, ...
           -eta .* sqrt_alpha_over_D .* sin_phi, ...
            xi .* sqrt_beta_over_D];
g.e_xi = [xi .* sqrt_beta_over_D .* cos_phi, ...
          xi .* sqrt_beta_over_D .* sin_phi, ...
          eta .* sqrt_alpha_over_D];
g.e_phi = [-sin_phi, cos_phi, zeros(size(phi))];
g.A_eta = A_eta;
g.A_eta_prime = A_eta_prime;
g.B = xi .* alpha;
g.B_prime = B_prime;
g.D_eta = D_eta;
g.D_xi = 2*xi;
g.nudged = nudged;
end
