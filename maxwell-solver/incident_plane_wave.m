function [E, H, khat, E0] = incident_plane_wave( ...
    k_m, n_m, theta_inc, phi_inc, pol, X, Y, Z)
%INCIDENT_PLANE_WAVE Plane wave for exp(-i*omega*t) convention.

validateattributes(k_m, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'k_m');
validateattributes(n_m, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'n_m');
validateattributes(theta_inc, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'theta_inc');
validateattributes(phi_inc, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'phi_inc');
validateattributes(pol, {'numeric'}, ...
    {'finite', 'vector', 'numel', 2}, mfilename, 'pol');
validateattributes(X, {'numeric'}, {'real', 'finite'}, mfilename, 'X');
validateattributes(Y, {'numeric'}, {'real', 'finite'}, mfilename, 'Y');
validateattributes(Z, {'numeric'}, {'real', 'finite'}, mfilename, 'Z');
if ~isequal(size(X), size(Y), size(Z))
    error('incident_plane_wave:InputSizeMismatch', ...
        'X, Y, and Z must have equal sizes.');
end

sin_theta = sin(theta_inc);
cos_theta = cos(theta_inc);
sin_phi = sin(phi_inc);
cos_phi = cos(phi_inc);
khat = [sin_theta*cos_phi, sin_theta*sin_phi, cos_theta];
e_te = [-sin_phi, cos_phi, 0];
e_tm = [cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta];
e1 = cos_phi*e_tm - sin_phi*e_te;
e2 = sin_phi*e_tm + cos_phi*e_te;
pol = pol(:).';
E0 = pol(1)*e1 + pol(2)*e2;
H0 = n_m * cross(khat, E0);

phase = exp(1i*k_m*(khat(1)*X(:) + khat(2)*Y(:) + khat(3)*Z(:)));
E = phase .* E0;
H = phase .* H0;
end
