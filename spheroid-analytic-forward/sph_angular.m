function W = sph_angular(eg, g)
%SPH_ANGULAR Evaluate regular metric-weighted angular combinations.

narginchk(2, 2);
required = {'mu', 'rmax', 'd'};
if ~isstruct(eg) || ~all(isfield(eg, required))
    error('sph_angular:InvalidEigenData', ...
        'eg must be a structure returned by sph_eigen.');
end
if ~isstruct(g) || ~isfield(g, 'eta') || ~isfield(g, 'beta')
    error('sph_angular:InvalidCoordinates', ...
        'g must contain eta and beta.');
end
eta = g.eta(:);
beta = g.beta(:);
if numel(eta) ~= numel(beta) || ~isreal(eta) || ~isreal(beta) || ...
        any(~isfinite(eta)) || any(~isfinite(beta)) || ...
        any(abs(eta) > 1) || any(beta < 0) || any(beta > 1)
    error('sph_angular:InvalidCoordinates', ...
        'eta and beta must be equally sized real values in [−1, 1] and [0, 1].');
end
mu = eg.mu;
r = (0:eg.rmax).';
if size(eg.d, 1) ~= numel(r)
    error('sph_angular:InvalidEigenData', ...
        'The row count of eg.d must equal eg.rmax + 1.');
end

ell_max = mu + eg.rmax;
U0 = reduced_basis(mu, ell_max, eta);
U1 = zeros(size(U0));
factor1 = zeros(size(r));
if eg.rmax >= 1
    U1(2:end, :) = reduced_basis(mu + 1, ell_max, eta);
    factor1(2:end) = sqrt(r(2:end) .* (2 * mu + r(2:end) + 1));
end
U2 = zeros(size(U0));
factor2 = zeros(size(r));
if eg.rmax >= 2
    U2(3:end, :) = reduced_basis(mu + 2, ell_max, eta);
    factor2(3:end) = sqrt(r(3:end) .* (2 * mu + r(3:end) + 1) .* ...
        (r(3:end) - 1) .* (2 * mu + r(3:end) + 2));
end

Stilde = eg.d.' * U0;
Stilde_prime = eg.d.' * (factor1 .* U1);
Stilde_second = eg.d.' * (factor2 .* U2);
eta = eta.';
beta = beta.';

if mu == 0
    W0 = Stilde;
    W1 = sqrt(beta) .* Stilde_prime;
    W2 = beta.^(3/2) .* Stilde_second;
    Wm = zeros(size(Stilde));
else
    base = beta.^((mu - 1) / 2);
    W0 = beta.^(mu / 2) .* Stilde;
    W1 = beta.^((mu + 1) / 2) .* Stilde_prime - ...
        mu * eta .* base .* Stilde;
    W2 = beta.^((mu + 3) / 2) .* Stilde_second - ...
        2 * mu * eta .* beta.^((mu + 1) / 2) .* Stilde_prime + ...
        (mu * (mu - 2) * eta.^2 - mu * beta) .* base .* Stilde;
    Wm = mu * base .* Stilde;
end

W = struct('W0', W0, 'W1', W1, 'W2', W2, 'Wm', Wm);
end

function U = reduced_basis(order, ell_max, eta)
count = ell_max - order + 1;
if count <= 0
    U = zeros(0, numel(eta));
    return
end
log_seed = 0.5 * (log((2 * order + 1) / 2) + ...
    gammaln(2 * order + 1) - 2 * order * log(2) - ...
    2 * gammaln(order + 1));
U = zeros(count, numel(eta));
U(1, :) = exp(log_seed);
if count == 1
    return
end
ell = order + 1;
a = sqrt((4 * ell^2 - 1) / (ell^2 - order^2));
U(2, :) = a * eta.' .* U(1, :);
for row = 3:count
    ell = order + row - 1;
    a = sqrt((4 * ell^2 - 1) / (ell^2 - order^2));
    b = sqrt(((ell - 1)^2 - order^2) / (4 * (ell - 1)^2 - 1));
    U(row, :) = a * (eta.' .* U(row - 1, :) - b * U(row - 2, :));
end
end
