function value = sph_semifocal(a, b)
%SPH_SEMIFOCAL Compute sqrt(abs(a^2-b^2)) without intermediate overflow.

validateattributes(a, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'a');
validateattributes(b, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'b');
a = double(a);
b = double(b);
scale = max(a, b);
value = scale*sqrt((abs(a - b)/scale)*(1 + min(a, b)/scale));
end
