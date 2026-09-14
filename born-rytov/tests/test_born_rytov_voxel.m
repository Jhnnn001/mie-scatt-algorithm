function test_born_rytov_voxel
% Regression: changing detector NA only filters the final complex fields.
forward_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(forward_dir);
N = 32;
dx = 0.1;
lambda = 1;
n_m = 1;
g = (-N/2:N/2-1)*dx;
[X,Y,Z] = ndgrid(g,g,g);
n = n_m + 0.12*exp(-((X-0.2).^2+2*Y.^2+Z.^2)/0.3^2);
theta = [0, asin(lambda/(N*dx))];
[b1,r1,u01,i1] = born_rytov_voxel( ...
    lambda,n,n_m,dx,theta,pi/2,0.4,2);
[b2,r2,u02,i2] = born_rytov_voxel( ...
    lambda,n,n_m,dx,theta,pi/2,0.8,2);
born_gap = 0;
rytov_gap = 0;
outside = 0;
for j = 1:numel(theta)
    B1 = fftshift(fft2(ifftshift(b1(:,:,j))));
    B2 = fftshift(fft2(ifftshift(b2(:,:,j))));
    R1 = fftshift(fft2(ifftshift(r1(:,:,j))));
    R2 = fftshift(fft2(ifftshift(r2(:,:,j))));
    born_gap = max(born_gap,norm(B1-B2.*i1.pupil,'fro')/norm(B1,'fro'));
    rytov_gap = max(rytov_gap,norm(R1-R2.*i1.pupil,'fro')/norm(R1,'fro'));
    outside = max(outside,norm(R1(~i1.pupil))/norm(R1,'fro'));
end
fprintf('NA regression: Born %.3e, Rytov %.3e, outside pupil %.3e\n', ...
    born_gap,rytov_gap,outside);
assert(born_gap < 1e-12,'Born detector pupils are inconsistent.');
assert(rytov_gap < 1e-12, ...
    'Changing detector NA changes the Rytov field before final filtering.');
assert(outside < 1e-12,'Rytov output contains frequencies outside detector NA.');
assert(isequal(u01,u02),'Detector NA changed the incident reference.');
assert(isequal(i1.rytov_phase,i2.rytov_phase), ...
    'The pre-objective Rytov phase must be independent of detector NA.');

[b0,r0] = born_rytov_voxel(lambda,ones(N,N,N),n_m,dx,theta,pi/2,0.4,2);
assert(all(b0 == 0,'all') && all(r0 == 0,'all'), ...
    'Zero contrast must give zero scattered fields.');
% A single voxel has a known point-quadrature Fourier amplitude and phase.
N = 16;
n = ones(N,N,N);
n(10,8,9) = 1.2;
[b,r,u0,info] = born_rytov_voxel(1,n,1,0.1,0,0,0.8,1);
[kx,ky] = ndgrid(info.kx,info.ky);
active = kx.^2+ky.^2 < (2*pi)^2;
kz = sqrt((2*pi)^2-kx(active).^2-ky(active).^2);
amplitude = (2*pi)^2*(1.2^2-1)*0.1^3;
S = 1i*amplitude./(2*kz).*exp(1i*kz) .* ...
    exp(-1i*(kx(active)*0.1-ky(active)*0.1));
[x,y] = ndgrid(info.x,info.y);
plane_waves = exp(1i*(x(:)*kx(active).'+y(:)*ky(active).'));
linear = reshape(plane_waves*S/(N*0.1)^2,N,N);
inside = kx(active).^2+ky(active).^2 <= (2*pi*0.8)^2;
expected_b = reshape(plane_waves*(S.*inside)/(N*0.1)^2,N,N);
expected_r = u0.*expm1(linear./u0);
expected_spectrum = fftshift(fft2(ifftshift(expected_r)));
expected_spectrum(~info.pupil) = 0;
expected_r = fftshift(ifft2(ifftshift(expected_spectrum)));
assert(norm(b-expected_b,'fro')/norm(expected_b,'fro') < 1e-12, ...
    'Point-source Born normalization or translation phase is incorrect.');
assert(norm(r-expected_r,'fro')/norm(expected_r,'fro') < 1e-12, ...
    'Point-source Rytov field is incorrect.');

expect_error(@() born_rytov_voxel(1,ones(8,8,8),1,0.6,0,0,0.1,3), ...
    'born_rytov_voxel:RytovNyquistViolation');
expect_error(@() born_rytov_voxel(1,ones(32,32,32),1,0.45,0.4,0,0.8,8), ...
    'born_rytov_voxel:InsufficientVoxelBandwidth');
fprintf('test_born_rytov_voxel: PASS\n');
end

function expect_error(action,identifier)
try
    action();
catch failure
    assert(strcmp(failure.identifier,identifier), ...
        'Expected %s, got %s.',identifier,failure.identifier);
    return
end
error('Expected error %s was not raised.',identifier);
end
