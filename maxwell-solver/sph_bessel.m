function [z, dz, d2z] = sph_bessel(kind, nu, x)
%SPH_BESSEL Spherical Bessel functions and first two derivatives.

kind = validate_kind(kind);
validateattributes(nu, {'numeric'}, ...
    {'real', 'finite', 'integer', 'nonnegative', 'vector'}, mfilename, 'nu');
validateattributes(x, {'numeric'}, {'finite', 'vector'}, mfilename, 'x');
nu = nu(:);
x = x(:).';
if ~strcmp(kind, 'j') && any(x == 0)
    error('sph_bessel:SingularAtZero', ...
        'Spherical y and h1 are singular at x=0.');
end

z = complex(zeros(numel(nu), numel(x)));
dz = z;
d2z = z;
nonzero = find(x ~= 0);
chunk_size = max(1, floor(1e6/max(1, numel(nu))));
for first = 1:chunk_size:numel(nonzero)
    columns = nonzero(first:min(first + chunk_size - 1, numel(nonzero)));
    xc = x(columns);
    [zc, znext] = values(kind, nu, xc);
    dzc = (nu./xc).*zc - znext;
    d2zc = (nu.*(nu - 1)./xc.^2 - 1).*zc + 2*znext./xc;
    z(:, columns) = zc;
    dz(:, columns) = dzc;
    d2z(:, columns) = d2zc;
end

zero = find(x == 0);
if ~isempty(zero)
    z(nu == 0, zero) = 1;
    dz(nu == 1, zero) = 1/3;
    d2z(nu == 0, zero) = -1/3;
    d2z(nu == 2, zero) = 2/15;
end
end

function [z, znext] = values(kind, nu, x)
[orders, arguments] = ndgrid(nu, x);
scale = sqrt(pi./(2*arguments));
switch kind
    case 'j'
        z = besselj(orders + 0.5, arguments).*scale;
        znext = besselj(orders + 1.5, arguments).*scale;
    case 'y'
        z = bessely(orders + 0.5, arguments).*scale;
        znext = bessely(orders + 1.5, arguments).*scale;
    otherwise
        z = besselh(orders + 0.5, 1, arguments).*scale;
        znext = besselh(orders + 1.5, 1, arguments).*scale;
end
end

function kind = validate_kind(kind)
if isstring(kind) && isscalar(kind)
    kind = char(kind);
end
if ~ischar(kind) || ~ismember(kind, {'j', 'y', 'h1'})
    error('sph_bessel:InvalidKind', 'kind must be j, y, or h1.');
end
end
