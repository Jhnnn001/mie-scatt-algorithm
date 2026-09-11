function test_gauss_legendre
% Focused checks for base-MATLAB Golub-Welsch quadrature.

[x, w] = gauss_legendre(4);
x_expected = [-0.8611363115940526; -0.3399810435848563; ...
               0.3399810435848563;  0.8611363115940526];
w_expected = [0.3478548451374539; 0.6521451548625461; ...
              0.6521451548625461; 0.3478548451374539];
assert(isequal(size(x), [4, 1]));
assert(isequal(size(w), [4, 1]));
assert(max(abs(x - x_expected)) < 5e-15);
assert(max(abs(w - w_expected)) < 5e-15);
assert(abs(sum(w .* x.^6) - 2/7) < 5e-15);

[x1, w1] = gauss_legendre(1);
assert(isequal(x1, 0) && isequal(w1, 2));

must_error(@() gauss_legendre(0));
must_error(@() gauss_legendre(2.5));
must_error(@() gauss_legendre([2, 3]));
must_error(@() gauss_legendre('4'));

fprintf('test_gauss_legendre: PASS\n');
end

function must_error(f)
did_error = false;
try
    f();
catch
    did_error = true;
end
assert(did_error);
end
