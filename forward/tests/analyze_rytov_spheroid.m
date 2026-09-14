function results = analyze_rytov_spheroid
%ANALYZE_RYTOV_SPHEROID Prolate baseline reported in the resume.
% Lengths in um; normal incidence, lab x polarization, detector NA=0.1.
here=fileparts(mfilename('fullpath'));
addpath(fileparts(here),fullfile(here,'..','..','spheroid-analytic-forward'));
lambda=.532; n_m=1.335381534; n_p=1.37;
a=5; b=2.5; N=128; dx=.1; padding=2; zd=7; NA=.1;
k0=2*pi/lambda; k=n_m*k0;
sol=spheroid_solve(lambda,n_p,n_m,a,b,0,0,[1,0], ...
    struct('tol',1e-6,'on_fail','error'));
assert(sol.info.validated,'Exact spheroidal solution failed validation.');
freq=(-N/2:N/2-1)*(2*pi/(N*dx)); [KX,KY]=ndgrid(freq,freq);
propagating=KX.^2+KY.^2<k^2; pupil=KX.^2+KY.^2<=(k0*NA)^2;
KZ=nan(N); KZ(propagating)=sqrt(k^2-KX(propagating).^2-KY(propagating).^2);
[exact,far_gap]=exact_spectrum(sol,KX,KY,KZ,pupil,zd);
xy=(-N/2:N/2-1)*dx;
inside=(xy(:).^2+xy.^2)/b^2+reshape(xy.^2,1,1,[])/a^2<=1;
n=n_m*ones(N,N,N); n(inside)=n_p;
[~,rytov,~,info]=born_rytov_voxel(lambda,n,n_m,dx,0,0,NA,zd, ...
    struct('pad_factor',padding));
S=fftshift(fft2(ifftshift(rytov)))*dx^2;
assert(norm(S(~pupil))/norm(S(:))<1e-12,'Rytov pupil leakage.');
born_error=relative(info.born_spectrum,exact);
corrected_rytov_error=relative(S,exact);
results=table(lambda,n_m,n_p,a,b,N,dx,padding,zd,NA,born_error, ...
    corrected_rytov_error,far_gap,'VariableNames', ...
    {'lambda_um','n_m','n_p','a_um','b_um','Nxy','dx_um','padding', ...
    'z_det_um','NA','born_error','corrected_rytov_error','exact_far_gap'});
writetable(results,fullfile(here,'rytov_analysis.csv'));
fprintf('PROLATE_BASELINE Born=%.6f%% Rytov=%.6f%%\n', ...
    100*born_error,100*corrected_rytov_error);
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
