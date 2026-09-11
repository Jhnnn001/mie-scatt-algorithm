function test_inverse
% Catches Fourier scaling/sign, phase handling, NA geometry and TV update bugs.
assert(exist('odt_prepare','file') == 2, 'Missing ODT data preparation.');
assert(exist('odt_reconstruct','file') == 2, 'Missing ODT reconstruction.');
rng(12);
lambda = 0.532; nm = 1.335381534; dx = 0.1;
n = nm*ones(16,16,16);
chi0 = 0.02;
n(9,9,9) = nm*sqrt(1+chi0);
theta = [0,0.12]; phi = [0,0.7];
[b,r,u0,g] = born_rytov_voxel(lambda,n,nm,dx,theta,phi,0.7,2);
db = odt_prepare(u0+b,u0,g,'born');
dr = odt_prepare(u0+r,u0,g,'rytov');
% A central voxel's continuous Fourier integral has no phase dependence.
assert(max(abs(db.samples-chi0*dx^3)) < 2e-12);
% Post-pupil logarithms are not identical to pre-pupil linear Rytov data.
rytov_gap=norm(dr.samples-db.samples)/norm(db.samples);
assert(rytov_gap>1e-9 && rytov_gap<.03);
scaled = odt_prepare((u0+b)*(2+3i),u0*(2+3i),g,'born');
assert(norm(scaled.samples-db.samples) < 1e-11);
% Analytic off-center point source: independent spectra, no voxel interpolation.
position = [0.4,-0.2,0.6];
[kx,ky] = ndgrid(g.kx,g.ky); p = g.pupil;
kz = sqrt(g.k^2-kx(p).^2-ky(p).^2);
exact_b = complex(zeros(size(b))); expected_q = []; expected_samples = [];
for j = 1:numel(theta)
    incident = g.k*[sin(theta(j))*cos(phi(j)),sin(theta(j))*sin(phi(j)),cos(theta(j))];
    q = [kx(p)-incident(1),ky(p)-incident(2),kz-incident(3)];
    sample = chi0*dx^3*exp(-1i*q*position.');
    spectrum = complex(zeros(size(p)));
    spectrum(p) = 1i*exp(1i*kz*g.z_det)./(2*kz)*g.k^2.*sample;
    exact_b(:,:,j) = fftshift(ifft2(ifftshift(spectrum)))/dx^2;
    expected_q = [expected_q;q]; %#ok<AGROW>
    expected_samples = [expected_samples;sample]; %#ok<AGROW>
end
shifted = odt_prepare(u0+exact_b,u0,g,'born');
assert(norm(shifted.q-expected_q,'fro') < 1e-12);
assert(norm(shifted.samples-expected_samples) < 1e-11);
% Wrapped phase exceeding pi, with a known zero-phase background border.
shape = sin(pi*(0:15)/15).^2;
phase = repmat(4.5*(shape.'*shape),1,1,2);
wrapped_U = u0.*exp(1i*phase);
automatic = odt_prepare(wrapped_U,u0,g,'rytov');
supplied = odt_prepare(wrapped_U,u0,g,'rytov',struct('unwrapped_phase',phase));
assert(norm(automatic.samples-supplied.samples) < 1e-11);
expect_error(@() odt_prepare(wrapped_U,u0,g,'rytov', ...
    struct('unwrapped_phase',phase+0.2)),'odt_prepare:InconsistentPhase');
expect_error(@() odt_prepare(u0*exp(2i),u0,g,'rytov'),'odt_prepare:AmbiguousPhase');
[b2,~,u02,g2] = born_rytov_voxel(lambda,n,nm,dx,theta,phi,0.4,2);
d2 = odt_prepare(u02+b2,u02,g2,'born');
assert(size(d2.q,1) < size(db.q,1));
[low_b,~,low_ref,low_g] = born_rytov_voxel(lambda,n,nm,dx,0,0,0.03,2);
low = odt_prepare(low_ref+low_b,low_ref,low_g,'born');
assert(all(low.q==0,'all'));
expect_error(@() odt_prepare(u0+b,zeros(size(u0)),g,'born'), ...
    'odt_prepare:InvalidReference');
expect_error(@() odt_prepare(zeros(size(u0)),u0,g,'rytov'), ...
    'odt_prepare:PhaseSingularity');
fprintf('Ewald scaling, Born/Rytov, reference and NA checks: PASS\n');

% Independent detector and volume dimensions must preserve physical samples.
small_geom=g; small_geom.grid_size=[8,8,8];
independent=odt_prepare(u0+b,u0,small_geom,'born');
assert(norm(independent.samples-db.samples)<1e-12);
assert(isequal(independent.operator.grid_size,[8,8,8]));
assert(exist('demo_inverse','file') == 2, 'Missing water/cell demonstration.');
demo = demo_inverse(0.55,'rytov',struct('grid_size',32,'a',1,'b',0.6, ...
    'n_angles',13,'plot',false,'max_iter',20,'outer_iter',2));
assert(numel(demo.reconstructions)==3 && all(isfinite(demo.metrics.relative_ri_error)));
assert(demo.n_m==nm && abs(demo.n_cell-1.33567771866)<1e-12 && demo.NA_det==0.55);
fprintf('inverse/tests/run_all: PASS\n');
end

function expect_error(action,id)
try
    action();
catch failure
    assert(strcmp(failure.identifier,id), ...
        'Expected %s, got %s.', id, failure.identifier);
    return
end
error('Expected error %s was not raised.',id);
end
