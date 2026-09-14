function [x, w] = gauss_legendre(Q)
%GAUSS_LEGENDRE Nodes and weights on [-1, 1] by Golub-Welsch.

validateattributes(Q, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'}, mfilename, 'Q');

if Q == 1
    x = 0;
    w = 2;
    return
end

j = (1:Q-1)';
b = j ./ sqrt(4*j.^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, order] = sort(diag(D));
w = 2 * V(1, order)'.^2;
end
