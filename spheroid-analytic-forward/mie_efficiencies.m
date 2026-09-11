function [Qext, Qsca, Qback, coeff] = mie_efficiencies(x, n_rel)
%MIE_EFFICIENCIES Lorenz-Mie efficiencies for exp(-i*omega*t).

validateattributes(x, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'x');
validateattributes(n_rel, {'numeric'}, {'finite', 'scalar'}, ...
    mfilename, 'n_rel');
if imag(n_rel) < 0 || n_rel == 0
    error('mie_efficiencies:InvalidRefractiveIndex', ...
        'n_rel must be nonzero with nonnegative imaginary part.');
end

n_size = max(x, abs(n_rel*x));
n0 = max(1, ceil(n_size + 4*n_size^(1/3) + 2));
converged = false;
for attempt = 1:8
    nmax = n0 + 10;
    [a, b, c, d] = coefficients(x, n_rel, nmax);
    q0 = efficiencies(x, a(1:n0), b(1:n0));
    q1 = efficiencies(x, a, b);
    finite = all(isfinite([a; b; c; d])) && all(isfinite([q0, q1]));
    if finite && max(abs(q1 - q0)./max(1, abs(q1))) < 1e-10
        converged = true;
        break
    end
    if ~finite
        break
    end
    n0 = nmax;
end
if ~converged
    error('mie_efficiencies:NoConvergence', ...
        'Mie coefficients or sums did not converge to finite values.');
end

Qext = q1(1);
Qsca = q1(2);
Qback = q1(3);
if nargout > 3
    coeff = struct('n', (1:nmax).', 'a', a, 'b', b, 'c', c, 'd', d);
end
end

function [a, b, c, d] = coefficients(x, m, nmax)
if m == 1
    a = zeros(nmax, 1);
    b = a;
    c = ones(nmax, 1);
    d = c;
    return
end

n = (1:nmax).';
psi_x = riccati_j(0:nmax, x);
psi_y = riccati_j(0:nmax, m*x);
xi_x = psi_x + 1i*riccati_y(0:nmax, x);
psi_x_prime = psi_x(1:end-1) - n.*psi_x(2:end)/x;
psi_y_prime = psi_y(1:end-1) - n.*psi_y(2:end)/(m*x);
xi_x_prime = xi_x(1:end-1) - n.*xi_x(2:end)/x;
psi_x = psi_x(2:end);
psi_y = psi_y(2:end);
xi_x = xi_x(2:end);

den_a = m*psi_y.*xi_x_prime - xi_x.*psi_y_prime;
den_b = psi_y.*xi_x_prime - m*xi_x.*psi_y_prime;
a = (m*psi_y.*psi_x_prime - psi_x.*psi_y_prime)./den_a;
b = (psi_y.*psi_x_prime - m*psi_x.*psi_y_prime)./den_b;
c = 1i*m./den_b;
d = 1i*m./den_a;
end

function psi = riccati_j(n, z)
psi = sqrt(pi*z/2)*besselj(n(:) + 0.5, z);
end

function chi = riccati_y(n, z)
chi = sqrt(pi*z/2)*bessely(n(:) + 0.5, z);
end

function q = efficiencies(x, a, b)
n = (1:numel(a)).';
w = 2*n + 1;
q = [2*sum(w.*real(a + b))/x^2, ...
    2*sum(w.*(abs(a).^2 + abs(b).^2))/x^2, ...
    abs(sum(w.*(-1).^n.*(a - b)))^2/x^2];
end
