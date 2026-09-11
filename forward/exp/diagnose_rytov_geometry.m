function results = diagnose_rytov_geometry
% Same axial thickness, wavelength, index, NA, and detector; wider transverse size.
out=fileparts(mfilename('fullpath'));
addpath(fullfile(out,'..','..','spheroid-analytic-forward'));
lambda=.532; nm=1.335381534; np=1.37; k0=2*pi/lambda; k=k0*nm;
N=1024; dx=.05; zd=7; freq=(-N/2:N/2-1)*(2*pi/(N*dx));
[KX,KY]=ndgrid(freq,freq); kt2=KX.^2+KY.^2; prop=kt2<k^2;
kz=zeros(N); kz(prop)=sqrt(k^2-kt2(prop)); pupil=kt2<(k0*.1)^2;
rows=struct([]);
for j=1:2
    a=5;
    if j==1
        b=2.5; s=scalar_prolate(lambda,np,nm,a,b,112);
        A=scalar_far_amplitude(s,kz(pupil)/k);
    else
        b=5; A=sphere_amplitude(k,k0*np,a,kz(pupil)/k);
    end
    exact=zeros(N); exact(pupil)=2i*pi*A.*exp(1i*kz(pupil)*zd)./kz(pupil);
    q=sqrt(b^2*kt2(prop)+a^2*(kz(prop)-k).^2);
    form=1-q.^2/10+q.^4/280; ii=q>1e-3;
    form(ii)=3*(sin(q(ii))-q(ii).*cos(q(ii)))./q(ii).^3;
    B=zeros(N); B(prop)=1i*k0^2*(np^2-nm^2)*(4*pi*a*b^2/3) ...
        *form.*exp(1i*kz(prop)*zd)./(2*kz(prop));
    psi=fftshift(ifft2(ifftshift(B)))/dx^2/exp(1i*k*zd);
    R=fftshift(fft2(ifftshift(exp(1i*k*zd)*expm1(psi))))*dx^2.*pupil;
    rows(j).a_um=a; rows(j).b_um=b;
    rows(j).central_phase=k0*(np-nm)*2*a;
    rows(j).rytov1_error=norm(R(:)-exact(:))/norm(exact(:));
    fprintf('Geometry a=%.6g b=%.3g, central phase=%.6g, R1 error %.4f%%\n', ...
        a,b,rows(j).central_phase,100*rows(j).rytov1_error);
end
results=struct2table(rows); writetable(results,fullfile(out,'geometry_control.csv'));
end

function A=sphere_amplitude(k,kp,a,eta)
% Independent exact scalar sphere partial waves avoid singular focal coordinates.
l=(0:112)'; x=k*a; y=kp*a;
jx=sqrt(pi/(2*x))*besselj(l+.5,x); jy=sqrt(pi/(2*y))*besselj(l+.5,y);
hx=sqrt(pi/(2*x))*besselh(l+.5,1,x);
dj=l/x.*jx-sqrt(pi/(2*x))*besselj(l+1.5,x);
di=l/y.*jy-sqrt(pi/(2*y))*besselj(l+1.5,y);
dh=l/x.*hx-sqrt(pi/(2*x))*besselh(l+1.5,1,x);
b=-(k*dj.*jy-kp*jx.*di)./(k*dh.*jy-kp*hx.*di);
assert(norm(b(end-5:end))<1e-10*norm(b),'Scalar sphere series tail is too large.');
P=zeros(numel(l),numel(eta)); P(1,:)=1; P(2,:)=eta.';
for n=1:numel(l)-2
    P(n+2,:)=((2*n+1)*eta.'.*P(n+1,:)-n*P(n,:))/(n+1);
end
A=(-1i/k*((2*l+1).*b).'*P).';
end
