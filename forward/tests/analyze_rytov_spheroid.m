function results = analyze_rytov_spheroid(case_indices)
% Exact vector Ex comparisons and independent continuous-shape diagnostics.
% Run from MATLAB: addpath('forward/tests'); analyze_rytov_spheroid
% Optional row indices rerun selected controls and return their table only.
% All lengths are um; normal incidence, lab x polarization, detector NA=0.1.
% No legacy spheroidal Born/Rytov solver is used.
tests_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(tests_dir),fullfile(tests_dir,'..','..','spheroid-analytic-forward'));
n_m = 1.335381534;
% name, lambda, n_p, a, b, Nx, Nz, dx, padding, z_det, use voxel FFT
cases = { ...
    'baseline',       .532,1.37,5,2.5,128,128,.1,2,7,true; ...
    'padding_4',      .532,1.37,5,2.5,128,128,.1,4,7,true; ...
    'dx_0p05',        .532,1.37,5,2.5,256,256,.05,2,7,true; ...
    'FOV_25p6',       .532,1.37,5,2.5,256,128,.1,2,7,true; ...
    'continuous_51',  .532,1.37,5,2.5,1024,0,.05,0,7,false; ...
    'continuous_102', .532,1.37,5,2.5,2048,0,.05,0,7,false; ...
    'near_surface',   .532,1.37,5,2.5,1024,0,.05,0,5.2,false; ...
    'farther_plane',  .532,1.37,5,2.5,1024,0,.05,0,14,false; ...
    'weak_contrast',  .532,n_m*1.001,5,2.5,128,128,.1,2,7,true; ...
    'half_size',      .532,1.37,2.5,1.25,128,128,.1,2,7,true; ...
    'lambda_633',     .633,1.37,5,2.5,128,128,.1,2,7,true; ...
    'half_scaled_z',  .532,1.37,2.5,1.25,1024,0,.05,0,3.5,false; ...
    'half_near',      .532,1.37,2.5,1.25,1024,0,.05,0,2.7,false; ...
    'lambda_633_fine',.633,1.37,5,2.5,1024,0,.05,0,7,false; ...
    'weak_fine',      .532,n_m*1.001,5,2.5,1024,0,.05,0,7,false};
write_results = nargin == 0;
if write_results, case_indices = 1:size(cases,1); end
cases = cases(case_indices,:);
last_key = [];
rows = struct([]);
for j = 1:size(cases,1)
    started = tic;
    [name,lambda,n_p,a,b,N,Nz,dx,padding,zd,use_voxel] = cases{j,:};
    k0 = 2*pi/lambda;
    k = n_m*k0;
    NA = .1;
    key = [lambda,n_p,a,b];
    if ~isequal(key,last_key)
        sol = spheroid_solve(lambda,n_p,n_m,a,b,0,0,[1,0], ...
            struct('tol',1e-6,'on_fail','error'));
        assert(sol.info.validated,'Exact spheroidal solution failed validation.');
        last_key = key;
    end
    freq = (-N/2:N/2-1)*(2*pi/(N*dx));
    [KX,KY] = ndgrid(freq,freq);
    propagating = KX.^2+KY.^2 < k^2;
    pupil = KX.^2+KY.^2 <= (k0*NA)^2;
    KZ = nan(N);
    KZ(propagating) = sqrt(k^2-KX(propagating).^2-KY(propagating).^2);
    [exact,far_gap] = exact_spectrum(sol,KX,KY,KZ,pupil,zd);
    % Closed-form FT of a homogeneous spheroid isolates voxel/interpolation error.
    Q = sqrt(b^2*(KX(propagating).^2+KY(propagating).^2) + ...
        a^2*(KZ(propagating)-k).^2);
    form = 1-Q.^2/10+Q.^4/280;
    regular = Q > 1e-3;
    qr = Q(regular);
    form(regular) = 3*(sin(qr)-qr.*cos(qr))./qr.^3;
    C = k0^2*(n_p^2-n_m^2)*(4*pi*a*b^2/3);
    linear_spectrum = complex(zeros(N));
    linear_spectrum(propagating) = 1i*C*form.*exp(1i*KZ(propagating)*zd) ./ ...
        (2*KZ(propagating));
    linear = fftshift(ifft2(ifftshift(linear_spectrum)))/dx^2;
    reference = exp(1i*k*zd);
    continuous_r = reference*expm1(linear/reference);
    continuous_spectrum = fftshift(fft2(ifftshift(continuous_r)))*dx^2;
    continuous_spectrum(~pupil) = 0;
    % A physical outgoing field must propagate linearly through homogeneous water.
    % Compare direct first Rytov at zd with Rytov formed near the exit then propagated.
    near_z = a+0.2;
    near_linear = linear_spectrum;
    near_linear(propagating) = near_linear(propagating) .* ...
        exp(1i*KZ(propagating)*(near_z-zd));
    near_linear = fftshift(ifft2(ifftshift(near_linear)))/dx^2;
    near_reference = exp(1i*k*near_z);
    near_r = near_reference*expm1(near_linear/near_reference);
    propagated = fftshift(fft2(ifftshift(near_r)))*dx^2;
    propagated(~pupil) = 0;
    propagated(pupil) = propagated(pupil).*exp(1i*KZ(pupil)*(zd-near_z));
    err_continuous = relative(continuous_spectrum,exact);
    err_born = relative(linear_spectrum.*pupil,exact);
    err_new = NaN; err_old = NaN; err_old_repupil = NaN;
    numerical_gap = NaN; volume_error = NaN; voxel_seconds = NaN;
    phi_peak = max(abs(linear/reference),[],'all');
    if use_voxel
        xy = (-N/2:N/2-1)*dx;
        zz = (-Nz/2:Nz/2-1)*dx;
        % Implicit expansion avoids three redundant 3-D coordinate volumes.
        inside = (xy(:).^2+xy.^2)/b^2 + reshape(zz.^2,1,1,[])/a^2 <= 1;
        n = n_m*ones(N,N,Nz);
        n(inside) = n_p;
        volume_error = nnz(inside)*dx^3/(4*pi*a*b^2/3)-1;
        [born,rytov,u0,info] = born_rytov_voxel( ...
            lambda,n,n_m,dx,0,0,NA,zd,struct('pad_factor',padding));
        clear n inside
        S = fftshift(fft2(ifftshift(rytov)))*dx^2;
        old = u0.*expm1(born./u0);
        old_spectrum = fftshift(fft2(ifftshift(old)))*dx^2;
        err_born = relative(info.born_spectrum,exact);
        err_new = relative(S,exact);
        err_old = relative(old_spectrum,exact);
        err_old_repupil = relative(old_spectrum.*pupil,exact);
        numerical_gap = relative(S,continuous_spectrum);
        phi_peak = max(abs(info.rytov_phase),[],'all');
        voxel_seconds = info.time;
        assert(norm(S(~pupil))/norm(S(:)) < 1e-12, ...
            'Corrected Rytov field leaked outside the pupil.');
        if strcmp(name,'weak_contrast')
            assert(err_new < .03,'Weak spheroid Rytov error exceeds 3%%.');
        end
    end
    rows(j).case_name = string(name);
    rows(j).lambda_um = lambda;
    rows(j).n_p = n_p;
    rows(j).a_um = a;
    rows(j).b_um = b;
    rows(j).Nxy = N;
    rows(j).Nz = Nz;
    rows(j).dx_um = dx;
    rows(j).padding = padding;
    rows(j).z_det_um = zd;
    rows(j).NA = NA;
    rows(j).pupil_modes = nnz(pupil);
    rows(j).phase_delay_estimate = k0*(n_p-n_m)*2*a;
    rows(j).born_error = err_born;
    rows(j).old_rytov_error = err_old;
    rows(j).old_repupil_error = err_old_repupil;
    rows(j).corrected_rytov_error = err_new;
    rows(j).continuous_rytov_error = err_continuous;
    rows(j).exit_propagated_error = relative(propagated,exact);
    rows(j).free_propagation_gap = relative(continuous_spectrum,propagated);
    rows(j).voxel_continuous_gap = numerical_gap;
    rows(j).voxel_volume_error = volume_error;
    rows(j).max_complex_phase = phi_peak;
    rows(j).exact_far_gap = far_gap;
    rows(j).voxel_seconds = voxel_seconds;
    rows(j).total_seconds = toc(started);
    fprintf('%-18s B=%7.3f%% oldR=%7.3f%% newR=%7.3f%% continuousR=%7.3f%% numerical=%6.3f%% (%.2fs)\n', ...
        name,100*err_born,100*err_old,100*err_new,100*err_continuous, ...
        100*numerical_gap,rows(j).total_seconds);
    results = struct2table(rows);
    if write_results
        writetable(results,fullfile(tests_dir,'rytov_analysis.csv'));
    end
end
end

function [spectrum,gap] = exact_spectrum(sol,KX,KY,KZ,pupil,zd)
% Outgoing far amplitude A relates to plane spectrum by U=2*pi*i*A/kz.
dirs = [KX(pupil),KY(pupil),KZ(pupil)]/sol.km;
radius = 1e6*max(sol.a,sol.b);
points = [radius*dirs;2*radius*dirs];
E = spheroid_eval(sol,points(:,1),points(:,2),points(:,3),'scattered');
m = nnz(pupil);
A1 = radius*exp(-1i*sol.km*radius)*E(1:m,1);
A2 = 2*radius*exp(-2i*sol.km*radius)*E(m+1:end,1);
gap = relative(A1,A2);
assert(gap < 1e-4,'Exact far-field extraction did not converge: %.3g.',gap);
spectrum = complex(zeros(size(pupil)));
spectrum(pupil) = 2*pi*1i./KZ(pupil).*A2.*exp(1i*KZ(pupil)*zd);
end

function e = relative(actual,expected)
e = norm(actual(:)-expected(:))/norm(expected(:));
end
