function results = diagnose_rytov(N)
% Isolate first-Rytov truncation from scalar-vector and voxel errors.
% Run: addpath('forward/exp'); diagnose_rytov
if nargin==0, N=512; end
out = fileparts(mfilename('fullpath'));
addpath(fullfile(out,'..','..','spheroid-analytic-forward'));
timer = tic;
lambda=.532; nm=1.335381534; np=1.37; a=5; b=2.5; NA=.1;
k0=2*pi/lambda; k=k0*nm; f0=k0^2*(np^2-nm^2);
dx=51.2/N; freq=(-N/2:N/2-1)*(2*pi/(N*dx));
[KX,KY]=ndgrid(freq,freq); kt2=KX.^2+KY.^2;
prop=kt2<k^2; pupil=kt2<=(k0*NA)^2;
KZ=zeros(N); KZ(prop)=sqrt(k^2-kt2(prop));
[eta,~,map]=unique(KZ(prop)/k);
Q=sqrt(b^2*kt2(prop)+a^2*(KZ(prop)-k).^2);
form=1-Q.^2/10+Q.^4/280; ii=Q>1e-3;
form(ii)=3*(sin(Q(ii))-Q(ii).*cos(Q(ii)))./Q(ii).^3;
born0=zeros(N); born0(prop)=1i*f0*(4*pi*a*b^2/3)*form./(2*KZ(prop));
scalar=scalar_prolate(lambda,np,nm,a,b,112);
refined=scalar_prolate(lambda,np,nm,a,b,132);
aa=scalar_far_amplitude(scalar,eta); ar=scalar_far_amplitude(refined,eta);
refinement_gap=rel(aa,ar); assert(refinement_gap<1e-7);
% Check the analytic far limit against finite-distance evaluation.
et=[.1;.6;.99]; radius=1e6*a;
pts=radius*[sqrt(1-et.^2),zeros(3,1),et];
[~,~,us]=scalar_prolate_eval(scalar,pts);
far_gap=rel(radius*exp(-1i*k*radius)*us,scalar_far_amplitude(scalar,et));
assert(far_gap<1e-3);
fprintf('Scalar L refinement %.3g; finite-radius check %.3g.\n',refinement_gap,far_gap);
vector=spheroid_solve(lambda,np,nm,a,b,0,0,[1,0],struct('tol',1e-6));
assert(vector.info.validated);
pts=2*radius*[KX(pupil),KY(pupil),KZ(pupil)]/k;
E=spheroid_eval(vector,pts(:,1),pts(:,2),pts(:,3),'scattered');
av=2*radius*exp(-2i*k*radius)*E(:,1);
vector0=zeros(N); vector0(pupil)=2i*pi*av./KZ(pupil);
base0=zeros(N); base0(prop)=2i*pi*aa(map)./KZ(prop);
vector_scalar_gap=rel(base0.*pupil,vector0);
fprintf('Exact scalar vs vector Ex on NA=.1: %.6f%%.\n',100*vector_scalar_gap);
% t scales scattering potential, not n-nm; this makes the Born coefficient linear.
% Symmetric small-contrast solves extract the first two exact Rytov coefficients.
h=.01; exact_coeff=cell(2,2);
for ih=1:2
    step=h/2^(ih-1); pm=cell(1,2);
    for sign_index=1:2
        t=(3-2*sign_index)*step;
        s=scalar_prolate(lambda,sqrt(nm^2+t*(np^2-nm^2)),nm,a,b,112);
        amp=scalar_far_amplitude(s,eta);
        pm{sign_index}=zeros(N); pm{sign_index}(prop)=2i*pi*amp(map)./KZ(prop);
    end
    exact_coeff{ih,1}=(pm{1}-pm{2})/(2*step);
    exact_coeff{ih,2}=(pm{1}+pm{2})/(2*step^2);
end
first_gap=rel(exact_coeff{2,1},born0);
second_gap=rel(exact_coeff{1,2},exact_coeff{2,2});
assert(first_gap<1e-3 && second_gap<1e-3);
metrics=table(refinement_gap,far_gap,vector_scalar_gap,first_gap,second_gap);
writetable(metrics,fullfile(out,sprintf('validation_N%d.csv',N)));
fprintf('Weak exact derivative vs Born %.3g; second coefficient h-refinement %.3g.\n',first_gap,second_gap);
rows=struct([]); j=0;
for t=[.1,.25,.5,1]
    if t==1
        exact0=base0;
    else
        s=scalar_prolate(lambda,sqrt(nm^2+t*(np^2-nm^2)),nm,a,b,112);
        amp=scalar_far_amplitude(s,eta);
        exact0=zeros(N); exact0(prop)=2i*pi*amp(map)./KZ(prop);
    end
    for zd=[5.2,7,14]
        phase=exp(1i*KZ*zd); u0=exp(1i*k*zd);
        phi1=ift(born0.*phase,dx)/u0;
        % psi2 = U2/U0 - (U1/U0)^2/2; use analytic Born for U1.
        phi2=ift(exact_coeff{2,2}.*phase,dx)/u0-.5*phi1.^2;
        R1=ft(u0*expm1(t*phi1),dx).*pupil;
        R2=ft(u0*expm1(t*phi1+t^2*phi2),dx).*pupil;
        exact=exact0.*phase.*pupil;
        j=j+1; rows(j).potential_scale=t; rows(j).z_um=zd;
        rows(j).born_error=rel(t*born0.*phase.*pupil,exact);
        rows(j).rytov1_error=rel(R1,exact);
        rows(j).rytov2_error=rel(R2,exact);
        rows(j).psi2_over_psi1=norm(t^2*phi2(:))/norm(t*phi1(:));
        rows(j).max_abs_psi1=max(abs(t*phi1),[],'all');
        rows(j).rytov1_vs_vector=NaN;
        if t==1, rows(j).rytov1_vs_vector=rel(R1,vector0.*phase); end
        fprintf('t=%.2f z=%4.1f B=%7.3f%% R1=%7.3f%% R2=%7.3f%% psi2/psi1=%.3g\n', ...
            t,zd,100*rows(j).born_error,100*rows(j).rytov1_error, ...
            100*rows(j).rytov2_error,rows(j).psi2_over_psi1);
    end
end
results=struct2table(rows);
writetable(results,fullfile(out,sprintf('contrast_distance_N%d.csv',N)));
save(fullfile(out,sprintf('spectral_N%d.mat',N)), ...
    'results','refinement_gap','far_gap','vector_scalar_gap', ...
    'first_gap','second_gap','phi1','phi2','dx','KZ','pupil','base0','born0');
fprintf('Spectral experiment completed in %.2fs.\n',toc(timer));
end

function u=ift(S,dx)
u=fftshift(ifft2(ifftshift(S)))/dx^2;
end

function S=ft(u,dx)
S=fftshift(fft2(ifftshift(u)))*dx^2;
end

function e=rel(a,b)
e=norm(a(:)-b(:))/norm(b(:));
end
