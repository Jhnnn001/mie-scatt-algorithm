function [T, dT, d2T] = sph_bessel_T(mu, r, x)
%SPH_BESSEL_T Evaluate j_(mu+r)(x)/x^mu and two derivatives.

validateattributes(mu, {'numeric'}, ...
    {'real', 'finite', 'integer', 'nonnegative', 'scalar'}, mfilename, 'mu');
validateattributes(r, {'numeric'}, ...
    {'real', 'finite', 'integer', 'nonnegative', 'vector'}, mfilename, 'r');
validateattributes(x, {'numeric'}, {'finite', 'vector'}, mfilename, 'x');
r = r(:);
x = x(:).';
tiny = abs(x) < 1e-3;
x_safe = x;
x_safe(tiny) = 1;

all_r = [r; r + 1; r + 2];
j = sph_bessel('j', mu + all_r, x_safe);
scaled = exp(log(j) - mu*log(x_safe));
count = numel(r);
T = scaled(1:count, :);
T1 = scaled(count + (1:count), :);
T2 = scaled(2*count + (1:count), :);
R = r;
X = x_safe;
dT = (R./X).*T - T1;
d2T = (R.*(R - 1)./X.^2).*T - ((2*R + 1)./X).*T1 + T2;

bad = repmat(tiny, count, 1) | invalid(j(1:count, :)) | ...
    invalid(j(count + (1:count), :)) | ...
    invalid(j(2*count + (1:count), :)) | ...
    invalid(T) | invalid(T1) | invalid(T2) | invalid(dT) | invalid(d2T);
[rows, columns] = find(bad);
for k = 1:numel(rows)
    [T(rows(k), columns(k)), dT(rows(k), columns(k)), ...
        d2T(rows(k), columns(k))] = power_series(mu, r(rows(k)), ...
        x(columns(k)));
end
end

function bad = invalid(value)
bad = value == 0 | ~isfinite(real(value)) | ~isfinite(imag(value));
end

function [value, derivative, second] = power_series(mu, r, x)
n = mu + r;
coefficient = exp(0.5*log(pi) - (n + 1)*log(2) - gammaln(n + 1.5));
if x == 0
    value = double(r == 0)*coefficient;
    derivative = double(r == 1)*coefficient;
    second = double(r == 2)*2*coefficient;
    if r == 0
        second = -coefficient/(2*n + 3);
    end
    return
end

term = exp(log(coefficient) + r*log(x));
value = term;
derivative = r*term/x;
second = r*(r - 1)*term/x^2;
scales = [abs(value), abs(derivative), abs(second)];
consecutive = 0;
for k = 1:max(100, ceil(2*abs(x) + 50))
    term = term*(-x^2)/(2*k*(2*n + 2*k + 1));
    power = r + 2*k;
    derivative_term = power*term/x;
    second_term = power*(power - 1)*term/x^2;
    value = value + term;
    derivative = derivative + derivative_term;
    second = second + second_term;
    scales = max(scales, [abs(term), abs(derivative_term), abs(second_term)]);
    converged = [abs(term), abs(derivative_term), abs(second_term)] <= ...
        8*eps*max(scales, realmin);
    if all(converged)
        consecutive = consecutive + 1;
        if consecutive == 3
            return
        end
    else
        consecutive = 0;
    end
end
error('sph_bessel_T:SeriesNoConvergence', ...
    'Scaled spherical-Bessel power series did not converge.');
end
