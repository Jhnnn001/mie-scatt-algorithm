function measurement = odt_exact_data(N, detector_dx, indices)
%ODT_EXACT_DATA Validated low-contrast oblate -> physical +z detector fields.
% Reuses the completed forward study, with its exact original polarization.
% No voxel truth, shape constraint, or inverse operator generates these fields.
if nargin<1, N=512; end
if nargin<2, detector_dx=.1; end
if nargin<3, indices=1:81; end
validateattributes(N,{'numeric'},{'integer','finite','scalar','>=',16});
validateattributes(detector_dx,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(indices,{'numeric'},{'integer','finite','vector','>=',1,'<=',81});
assert(mod(N,2)==0 && numel(unique(indices))==numel(indices),'Even detector, unique directions required.');
root=fileparts(fileparts(mfilename('fullpath')));
reference_dir=fullfile(root,'experiments','contrast-sweep');
addpath(fullfile(root,'experiments','field-analysis'),fullfile(root,'maxwell-solver'));
checks=readtable(fullfile(reference_dir,'contrast_1_percent_reference_validation.csv'));
assert(height(checks)==1 && checks.case_id==3 && ...
    all(isfinite(checks{1,2:end})) && max(checks{1,2:5})<1e-6 && ...
    checks.public_far_gap<1e-5,'Unvalidated exact reference.');
loaded=load(fullfile(reference_dir,'contrast_1_percent_reference.mat'),'s'); s=loaded.s;
assert(s.a==3 && s.b==5 && s.NA_det==1 && s.z_det==5 && ...
    abs(s.n_m-1.335381534)<1e-12 && abs(s.n_p-1.33567771866)<1e-12, ...
    'The selected reference is not the approved delta-n/100 system.');
angles=readtable(fullfile(reference_dir,'contrast_1_percent_angles.csv'));
assert(height(angles)==81 && all(angles.case_id==3) && ...
    all(isfinite(angles{:,2:end}),'all') && ...
    height(unique(angles(:,{'illumination_NA','phi_deg'})))==81, ...
    'Incomplete or invalid illumination data.');
all_weights=zeros(81,1);
api=oblate_field_analysis('helpers'); all_weights(:)=api.weights(angles);
angles=angles(indices,:); weights=all_weights(indices); weights=weights/sum(weights);
k=s.km; k0=2*pi/s.lambda;
assert(pi/detector_dx>k0*s.NA_det,'Detector Nyquist violation.');
f=(-N/2:N/2-1)*2*pi/(N*detector_dx); [KX,KY]=ndgrid(f,f);
pupil=KX.^2+KY.^2<=(k0*s.NA_det)^2;
prop=KX.^2+KY.^2<k^2; KZ=nan(N); KZ(prop)=sqrt(k^2-KX(prop).^2-KY(prop).^2);
eta=KZ(pupil)/k; az=atan2(KY(pupil),KX(pupil));
factor=2*pi*1i*exp(1i*KZ(pupil)*s.z_det)./KZ(pupil);
x=(-N/2:N/2-1)*detector_dx; [X,Y]=ndgrid(x,x);
ni=height(angles); U=complex(zeros(N,N,ni)); u0=U; scalar_U=U;
co_spectrum=complex(zeros(nnz(pupil),ni)); scalar_spectrum=co_spectrum;
ki=zeros(ni,3); timer=tic;
for phi_deg=unique(angles.phi_deg).'
    phi=phi_deg*pi/180;
    selected=find(angles.phi_deg==phi_deg);
    theta_indices=round(angles.illumination_NA(selected)/.05)+1;
    sub=api.subset(s,theta_indices.');
    [EV,US]=api.far_fields(sub,eta,az-phi);
    for local=1:numel(selected)
        j=selected(local); theta=angles.theta_deg(j)*pi/180;
        e0=cos(phi)*[cos(theta),0,-sin(theta)]-sin(phi)*[0,1,0];
        exact=cos(phi)*EV(:,:,2*local-1)-sin(phi)*EV(:,:,2*local);
        co_spectrum(:,j)=(exact*e0.').*factor;
        scalar_spectrum(:,j)=US(:,local).*factor;
        ki(j,:)=k*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
        ref=exp(1i*(ki(j,1)*X+ki(j,2)*Y+ki(j,3)*s.z_det));
        spectrum=complex(zeros(N)); spectrum(pupil)=co_spectrum(:,j);
        U(:,:,j)=ref+fftshift(ifft2(ifftshift(spectrum)))/detector_dx^2;
        spectrum(pupil)=scalar_spectrum(:,j);
        scalar_U(:,:,j)=ref+fftshift(ifft2(ifftshift(spectrum)))/detector_dx^2;
        u0(:,:,j)=ref;
    end
    fprintf('Exact detector N=%d dx=%g phi=%g: %d directions, %.1fs\n', ...
        N,detector_dx,phi_deg,numel(selected),toc(timer));
end
geometry=struct('grid_size',[128,128,80],'dx',.1,'detector_dx',detector_dx, ...
    'k',k,'k0',k0,'kx',f,'ky',f,'kz',KZ,'pupil',pupil, ...
    'k_incident',ki,'z_det',s.z_det,'NA_det',s.NA_det,'pad_factor',2);
measurement=struct('U',U,'u0',u0,'scalar_U',scalar_U, ...
    'co_spectrum',co_spectrum,'scalar_spectrum',scalar_spectrum, ...
    'geometry',geometry,'weights',weights,'angles',angles,'parameters',s.parameters, ...
    'indices',indices,'reference_validation',checks,'elapsed',toc(timer));
end
