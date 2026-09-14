function [uB,uR,spectra] = odt_predict(chi, geom, padding)
%ODT_PREDICT Same voxel Born/Rytov equations on an independent detector grid.
% The small extension needed beyond born_rytov_voxel is detector/volume size
% independence. FFT scale, interpn and post-exponential pupil are unchanged.
if nargin<3, padding=2; end
validateattributes(chi,{'numeric'},{'real','finite','nonempty'});
sz=geom.grid_size; assert(isequal(size(chi),sz),'Wrong volume size.');
ps=sz*padding; dx=geom.dx; det_dx=dx;
if isfield(geom,'detector_dx'), det_dx=geom.detector_dx; end
if pi/det_dx<=geom.k
    error('odt_predict:Bandwidth','Resolve the pre-pupil hemisphere with detector_dx < lambda/(2*n_m).');
end
fx=(-ps(1)/2:ps(1)/2-1)*2*pi/(ps(1)*dx);
fy=(-ps(2)/2:ps(2)/2-1)*2*pi/(ps(2)*dx);
[KX,KY]=ndgrid(geom.kx,geom.ky); N=size(KX);
prop=KX.^2+KY.^2<geom.k^2; pupil=geom.pupil;
KZ=zeros(N); KZ(prop)=sqrt(geom.k^2-KX(prop).^2-KY(prop).^2);
qz=[min(KZ(prop))-max(geom.k_incident(:,3)),max(KZ(prop))-min(geom.k_incident(:,3))];
zc=qz/(2*pi/(ps(3)*dx))+ps(3)/2+1;
assert(all(zc>=1 & zc<=ps(3)),'Pre-pupil Ewald samples exceed volume z bandwidth.');
planes=max(1,min(ps(3)-1,floor(zc(1)))):min(ps(3),floor(zc(2))+1);
fz=(planes-ps(3)/2-1)*2*pi/(ps(3)*dx);
z=(-sz(3)/2:sz(3)/2-1)*dx;
axial=reshape(reshape(chi,[],sz(3))*exp(-1i*z(:)*fz),[sz(1:2),numel(planes)]);
padded=complex(zeros([ps(1:2),numel(planes)])); offset=(ps-sz)/2;
padded(offset(1)+(1:sz(1)),offset(2)+(1:sz(2)),:)=axial;
F=fftshift(fftshift(fft2(ifftshift(ifftshift(padded,1),2)),1),2)*(dx^3*geom.k^2);
clear padded axial
factor=1i*exp(1i*KZ(prop)*geom.z_det)./(2*KZ(prop));
[X,Y]=ndgrid((-N(1)/2:N(1)/2-1)*det_dx,(-N(2)/2:N(2)/2-1)*det_dx);
ni=size(geom.k_incident,1); uB=complex(zeros([N,ni])); uR=uB;
bs=complex(zeros(nnz(pupil),ni)); rs=bs;
for j=1:ni
    ki=geom.k_incident(j,:);
    query=[KX(prop)-ki(1),KY(prop)-ki(2),KZ(prop)-ki(3)];
    if any(query<[fx(1),fy(1),fz(1)] | query>[fx(end),fy(end),fz(end)],'all')
        error('odt_predict:Bandwidth','Pre-pupil Ewald samples exceed volume bandwidth.');
    end
    B=complex(zeros(N));
    B(prop)=factor.*interpn(fx,fy,fz,F,query(:,1),query(:,2),query(:,3),'linear');
    pre=fftshift(ifft2(ifftshift(B)))/det_dx^2;
    ref=exp(1i*(ki(1)*X+ki(2)*Y+ki(3)*geom.z_det));
    R=fftshift(fft2(ifftshift(ref.*expm1(pre./ref))))*det_dx^2;
    R(~pupil)=0; B(~pupil)=0;
    bs(:,j)=B(pupil); rs(:,j)=R(pupil);
    uB(:,:,j)=fftshift(ifft2(ifftshift(B)))/det_dx^2;
    uR(:,:,j)=fftshift(ifft2(ifftshift(R)))/det_dx^2;
end
spectra=struct('born',bs,'rytov',rs);
end
