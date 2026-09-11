%EXAMPLE_XZ_SLICE Total-field magnitude for Barton's prolate n = 1.33 case.

root_directory = fileparts(mfilename('fullpath'));
addpath(root_directory);

semifocal = 1;
xi0 = 1.154701;
a = xi0*semifocal;
b = sqrt(xi0^2 - 1)*semifocal;
h_ext = 27.494593;
lambda = 2*pi/h_ext;

sol = spheroid_solve(lambda, 1.33, 1, a, b, pi/6, 0, [0, 1]);
axis_values = linspace(-1.5, 1.5, 301)*semifocal;
[X, Z] = meshgrid(axis_values, axis_values);
Y = zeros(size(X));
E = spheroid_eval(sol, X, Y, Z, 'total');
field_magnitude = reshape(sqrt(sum(abs(E).^2, 2)), size(X));

if usejava('jvm')
    imagesc(axis_values/semifocal, axis_values/semifocal, field_magnitude);
    set(gca, 'YDir', 'normal');
    axis image;
    xlabel('x / l');
    ylabel('z / l');
    title('Barton prolate case, n_p = 1.33: |E|');
    colorbar;
else
    fprintf('example_xz_slice: field computed; plotting skipped without JVM.\n');
end
