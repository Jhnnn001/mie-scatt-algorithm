function api = oblate_na_sweep(stage)
%OBLATE_NA_SWEEP Reproducible Born/Rytov experiment, lengths in um.
% Run stages 'check', 'reference', 'sweep', 'fieldcheck', 'report' in that order.
% Experimental batch assembly uses the public solver's bases and equations;
% its validation is recorded separately from spheroid_solve/info.validated.
if nargin==0, stage='check'; end
api=[];
out=fileparts(mfilename('fullpath'));
addpath(fileparts(out),fullfile(out,'..','..','spheroid-analytic-forward'));
p=struct('lambda',.532,'n_m',1.335381534,'n_p',1.365,'a',3,'b',5, ...
    'NA_det',1,'z_det',5,'illumination_NA',(0:.05:.5).');
p.theta=asin(p.illumination_NA/p.n_m);
switch stage
    case 'helpers'
        % Reuse the same experimental equations in the size/contrast studies.
        api=struct('reference',@batch_reference,'far_fields',@far_fields, ...
            'public_far_check',@public_far_check,'subset',@subset_reference, ...
            'spectra',@continuous_spectra,'weights',@illumination_weights, ...
            'near_fields',@near_fields,'scalar_scaled_reference',@scalar_scaled_reference);
    case 'check'
        q=p; q.lambda=2*pi; q.n_m=1; q.n_p=1+1e-5;
        q.a=.6; q.b=1; q.theta=[0;.3]; q.NA_det=.8;
        s=batch_reference(q,18,12);
        eta=linspace(-.9,.99,21).'; az=(0:20).'*1.1;
        [E,U]=far_fields(s,eta,az);
        dirs=[sqrt(1-eta.^2).*cos(az),sqrt(1-eta.^2).*sin(az),eta];
        for j=1:numel(q.theta)
            ki=[sin(q.theta(j)),0,cos(q.theta(j))];
            Q=sqrt(q.b^2*sum((dirs(:,1:2)-ki(1:2)).^2,2)+ ...
                q.a^2*(dirs(:,3)-ki(3)).^2);
            B=(q.n_p^2-q.n_m^2)*(4*pi*q.a*q.b^2/3)*form_factor(Q)/(4*pi);
            assert(relative(U(:,j),B)<1e-4,'Scalar weak-contrast Born check.');
            for pol=1:2
                [~,~,~,e0]=incident_plane_wave(1,1,q.theta(j),0, ...
                    double([pol==1,pol==2]),0,0,0);
                BV=B.*(e0-(dirs*e0.').*dirs);
                assert(relative(E(:,:,2*j-2+pol),BV)<1e-4, ...
                    'Vector weak-contrast Born check.');
            end
        end
        gap=public_far_check(s);
        assert(gap<1e-5,'Public evaluator far-limit check.');
        [~,~,~,ex]=incident_plane_wave(1,1,.3,pi/2,[1,0],0,0,0);
        assert(norm(ex-[1,0,0])<1e-14,'Original polarization convention.');
        fprintf('OBLATE_IMPLEMENTATION_CHECKS_PASS public_far_gap=%.3g\n',gap);
    case 'reference'
        for fine=0:1
            filename=fullfile(out,sprintf('oblate_reference_%d.mat',fine));
            if isfile(filename) && ~isempty(whos('-file',filename,'s'))
                loaded=load(filename,'s'); s=loaded.s;
                assert(isequal(s.parameters,p),'Reference checkpoint parameters differ.');
                fprintf('Loaded %s\n',filename);
            else
                s=batch_reference(p,109+5*fine,55+4*fine);
                s.public_far_gap=public_far_check(s);
                assert(s.public_far_gap<1e-5,'Far limit vs public evaluator.');
                % Far-field checkpoints need no dense outgoing ODE trajectories.
                % Re-run this stage to obtain a full public-evaluator structure.
                for mu=1:numel(s.mu_data)
                    s.mu_data{mu}.rod_ext=rmfield(s.mu_data{mu}.rod_ext,'solution');
                end
                save(filename,'s','-v7');
            end
        end
        low=load(fullfile(out,'oblate_reference_0.mat'),'s');
        [q,w]=gauss_legendre(96); q=(q+1)/2; w=w/2;
        [Q,Phi]=ndgrid(q,2*pi*(0:127)/128);
        eta=sqrt(1-(Q/p.n_m).^2);
        [E0,U0]=far_fields(low.s,eta(:),Phi(:));
        [E1,U1]=far_fields(s,eta(:),Phi(:));
        weights=repmat(w.*q./(p.n_m^2-q.^2),128,1);
        cutoff_vector=relative(E0.*sqrt(weights),E1.*sqrt(weights));
        cutoff_scalar=relative(U0.*sqrt(weights),U1.*sqrt(weights));
        assert(max(cutoff_vector,cutoff_scalar)<1e-6,'Reference cutoff convergence.');
        validation=table(cutoff_vector,cutoff_scalar,s.vector_boundary_error, ...
            s.scalar_boundary_error,s.public_far_gap,'VariableNames', ...
            {'cutoff_vector','cutoff_scalar','vector_boundary','scalar_boundary','public_far_gap'});
        writetable(validation,fullfile(out,'oblate_reference_validation.csv'));
        fprintf('OBLATE_REFERENCES_PASS vector %.3g scalar %.3g\n',cutoff_vector,cutoff_scalar);
    case 'sweep'
        load(fullfile(out,'oblate_reference_1.mat'),'s');
        assert(isequal(s.parameters,p),'Reference parameters differ.');
        run_sweep(s,out);
    case 'report'
        summarize_results(out);
    case 'fieldcheck'
        field_convergence(p,out);
    case {'xz','xz_validate'}
        p.theta=[0;5;10;15]*pi/180;
        p.illumination_NA=p.n_m*sin(p.theta);
        s=batch_reference(p,114,59);
        if strcmp(stage,'xz'), run_xz(s,out); else, validate_xz(s,out); end
    otherwise
        error('Unknown experiment stage.');
end
end

function s=batch_reference(p,L,M)
% One matrix per signed m, multiple incidence/polarization right-hand sides.
s=p; s.parameters=p; s.L=L; s.M=M; s.type='oblate'; s.route='spheroid';
s.semifocal=sqrt(p.b^2-p.a^2); s.xi0=p.a/s.semifocal;
s.km=2*pi*p.n_m/p.lambda; s.c_ext=s.km*s.semifocal;
s.c_int=2*pi*p.n_p/p.lambda*s.semifocal;
s.Qtheta=ceil(L+max(s.c_ext,s.c_int)+30); s.Qphi=256;
[g,w]=meridian(s,s.Qtheta); [gv,wv]=meridian(s,s.Qtheta+17);
inc=boundary_incident(s,g); iv=boundary_incident(s,gv);
n=numel(p.theta); nv=2*n;
s.mu_data=cell(M+1,1);
template=struct('m',0,'mu',0,'a',[],'b',[],'c',[],'d',[], ...
    'so',[],'si',[],'analytically_zeroed',false);
s.modes=repmat(template,2*M+1,1);
errv=zeros(1,nv); errs=zeros(1,n);
scalev=repelem([1;1;1/p.n_m;1/p.n_m;1/p.n_m^2;1/p.n_m],numel(wv));
wv6=repmat(wv,6,1); wv2=repmat(wv,2,1);
timer=tic;
for mu=0:M
    eo=sph_eigen('oblate',mu,s.c_ext,L);
    ei=sph_eigen('oblate',mu,s.c_int,L);
    [rod,ok]=sph_radial_ode(eo,s.c_ext,s.xi0,1.5);
    assert(ok.ok,'Outgoing radial ODE check failed.');
    s.mu_data{mu+1}=struct('mu',mu,'eg_ext',eo,'eg_int',ei,'rod_ext',rod);
    Wo=sph_angular(eo,g); Wi=sph_angular(ei,g);
    Vo=sph_radial_ode_eval(rod,g); [Vi,ok]=sph_radial(ei,1,g,s.c_int);
    assert(ok,'Regular radial check failed.');
    Wov=sph_angular(eo,gv); Wiv=sph_angular(ei,gv);
    Vov=sph_radial_ode_eval(rod,gv); [Viv,ok]=sph_radial(ei,1,gv,s.c_int);
    assert(ok,'Regular radial validation check failed.');
    S=scalar_matrix(Wo,Vo,Wi,Vi,g,s.c_ext);
    Sv=scalar_matrix(Wov,Vov,Wiv,Viv,gv,s.c_ext);
    project=kron(eye(4),Wo.W0.*w.');
    ms=unique([-mu,mu]); count=L-mu+1;
    for m=ms
        bin=mod(m,s.Qphi)+1;
        C=vector_matrix(Wo,Vo,Wi,Vi,g,m,s);
        rhs=-reshape(inc.v(:,bin,:),[],nv);
        coefficients=scaled_solve(project*C(1:4*numel(w),:), ...
            project*rhs(1:4*numel(w),:));
        srhs=-reshape(inc.s(:,bin,:),[],n);
        sc=scaled_solve(S.*sqrt(repmat(w,2,1)),srhs.*sqrt(repmat(w,2,1)));
        k=m+M+1;
        s.modes(k)=struct('m',m,'mu',mu,'a',coefficients(1:count,:), ...
            'b',coefficients(count+(1:count),:),'c',coefficients(2*count+(1:count),:), ...
            'd',coefficients(3*count+(1:count),:),'so',sc(1:count,:), ...
            'si',sc(count+1:end,:),'analytically_zeroed',false);
        Cv=vector_matrix(Wov,Vov,Wiv,Viv,gv,m,s);
        rv=(Cv*coefficients+reshape(iv.v(:,bin,:),[],nv)).*scalev;
        rs=Sv*sc+reshape(iv.s(:,bin,:),[],n);
        errv=errv+sum(wv6.*abs(rv).^2,1);
        errs=errs+sum(wv2.*abs(rs).^2,1);
    end
    fprintf('Reference L=%d M=%d mu=%d elapsed=%.1fs\n',L,M,mu,toc(timer));
end
used=mod(-M:M,s.Qphi)+1; unused=setdiff(1:s.Qphi,used);
denv=squeeze(sum(sum(wv6.*abs(iv.v.*scalev).^2,1),2)).';
dens=squeeze(sum(sum(wv2.*abs(iv.s).^2,1),2)).';
errv=errv+squeeze(sum(sum(wv6.*abs(iv.v(:,unused,:).*scalev).^2,1),2)).';
errs=errs+squeeze(sum(sum(wv2.*abs(iv.s(:,unused,:)).^2,1),2)).';
s.vector_boundary_by_case=sqrt(errv./denv);
s.scalar_boundary_by_case=sqrt(errs./dens);
s.vector_boundary_error=max(s.vector_boundary_by_case);
s.scalar_boundary_error=max(s.scalar_boundary_by_case);
assert(s.vector_boundary_error<1e-6 && s.scalar_boundary_error<1e-6, ...
    'Independent boundary residual: vector %.3g scalar %.3g.', ...
    s.vector_boundary_error,s.scalar_boundary_error);
s.solve_seconds=toc(timer);
fprintf('BATCH_BOUNDARY_PASS vector=%.3g scalar=%.3g time=%.1fs\n', ...
    s.vector_boundary_error,s.scalar_boundary_error,s.solve_seconds);
end

function [g,w]=meridian(s,N)
[x,w]=gauss_legendre(N); t=pi*(x+1)/2; w=pi*w.*sin(t)/2;
g=sph_coords('oblate',s.b/s.semifocal*sin(t),zeros(size(t)),s.xi0*cos(t));
end

function inc=boundary_incident(s,g)
eta=g.eta; n=numel(s.theta); nq=numel(eta); phi=2*pi*(0:s.Qphi-1)/s.Qphi;
[Eta,Phi]=ndgrid(eta,phi); rho=s.b*sqrt(1-Eta.^2);
X=rho.*cos(Phi); Y=rho.*sin(Phi); Z=s.a*Eta;
gg=sph_coords('oblate',X(:)/s.semifocal,Y(:)/s.semifocal,Z(:)/s.semifocal);
inc.v=complex(zeros(6*nq,s.Qphi,2*n));
inc.s=complex(zeros(2*nq,s.Qphi,n));
for j=1:n
    for pol=1:2
        [E,H]=incident_plane_wave(s.km,s.n_m,s.theta(j),0,double([pol==1,pol==2]),X,Y,Z);
        fields={sum(E.*gg.e_eta,2),sum(E.*gg.e_phi,2), ...
            sum(H.*gg.e_eta,2),sum(H.*gg.e_phi,2), ...
            s.n_m^2*sum(E.*gg.e_xi,2),sum(H.*gg.e_xi,2)};
        for k=1:6
            inc.v((k-1)*nq+(1:nq),:,2*j-2+pol)=fft(reshape(fields{k},nq,[]),[],2)/s.Qphi;
        end
    end
    U=exp(1i*s.km*(X*sin(s.theta(j))+Z*cos(s.theta(j))));
    dU=1i*(s.xi0/sqrt(s.xi0^2+1)*sqrt(1-Eta.^2).*cos(Phi)*sin(s.theta(j)) ...
        +Eta*cos(s.theta(j))).*U;
    inc.s(:,:,j)=[fft(U,[],2);fft(dU,[],2)]/s.Qphi;
end
end

function C=vector_matrix(Wo,Vo,Wi,Vi,g,m,s)
[Mo,No]=sph_vecwave(Wo,Vo,g,m,s.c_ext);
[Mi,Ni]=sph_vecwave(Wi,Vi,g,m,s.c_int);
C=[No.eta.',Mo.eta.',-Ni.eta.',-Mi.eta.'; ...
   No.phi.',Mo.phi.',-Ni.phi.',-Mi.phi.'; ...
   -1i*s.n_m*Mo.eta.',-1i*s.n_m*No.eta.',1i*s.n_p*Mi.eta.',1i*s.n_p*Ni.eta.'; ...
   -1i*s.n_m*Mo.phi.',-1i*s.n_m*No.phi.',1i*s.n_p*Mi.phi.',1i*s.n_p*Ni.phi.'; ...
   s.n_m^2*No.xi.',s.n_m^2*Mo.xi.',-s.n_p^2*Ni.xi.',-s.n_p^2*Mi.xi.'; ...
   -1i*s.n_m*Mo.xi.',-1i*s.n_m*No.xi.',1i*s.n_p*Mi.xi.',1i*s.n_p*Ni.xi.'];
end

function C=scalar_matrix(Wo,Vo,Wi,Vi,g,co)
C=[(Wo.W0.*Vo.V0).',-(Wi.W0.*Vi.V0).'; ...
   (Wo.W0.*Vo.V1).'/sqrt(g.alpha(1))/co, ...
   -(Wi.W0.*Vi.V1).'/sqrt(g.alpha(1))/co];
end

function c=scaled_solve(A,b)
scales=max(abs(A),[],1); assert(all(scales>0 & isfinite(scales)));
A=A./scales;
if size(A,1)==size(A,2) && rcond(A)<=size(A,1)*eps
    c=lsqminnorm(A,b)./scales.';
else
    c=(A\b)./scales.';
end
assert(all(isfinite(c),'all'),'Nonfinite boundary coefficients.');
end

function [E,U]=far_fields(s,eta,phi)
% Analytic r*exp(-ikr) limit of the same waves used by spheroid_eval.
eta=eta(:); phi=phi(:); [eu,~,map]=unique(eta);
g=struct('eta',eu,'beta',max(0,1-eu.^2)); n=numel(s.theta);
Ae=complex(zeros(numel(eta),2*n)); Ap=Ae; U=complex(zeros(numel(eta),n));
do_vector=~isfield(s,'scalar_only') || ~s.scalar_only;
for mu=0:s.M
    d=s.mu_data{mu+1}; W=sph_angular(d.eg_ext,g); r=d.rod_ext;
    f=s.semifocal*exp(-1i*r.c*r.xi0)*exp(-r.qtilde_xi0(:)) .* ...
        r.asymptotic_coefficients(1,:).' ./ ...
        (r.asymptotic_argument_scale*r.asymptotic_U_start(:));
    for k=find([s.modes.mu]==mu)
        mode=s.modes(k); phase=exp(1i*mode.m*phi); sg=sign(mode.m);
        if do_vector
            ae=(1i*(mode.a.*f).'*W.W1-1i*sg*(mode.b.*f).'*W.Wm).';
            ap=(-sg*(mode.a.*f).'*W.Wm+(mode.b.*f).'*W.W1).';
            Ae=Ae+phase.*ae(map,:); Ap=Ap+phase.*ap(map,:);
        end
        us=((mode.so.*f).'*W.W0).';
        U=U+phase.*us(map,:);
    end
end
E=complex(zeros(numel(eta),3,2*n));
E(:,1,:)=reshape(-eta.*cos(phi).*Ae-sin(phi).*Ap,numel(eta),1,[]);
E(:,2,:)=reshape(-eta.*sin(phi).*Ae+cos(phi).*Ap,numel(eta),1,[]);
E(:,3,:)=reshape(sqrt(max(0,1-eta.^2)).*Ae,numel(eta),1,[]);
end

function gap=public_far_check(s)
eta=[.31;.61;.81;.93;.999]; az=[.1;.9;1.8;3.1;4.5];
dirs=[sqrt(1-eta.^2).*cos(az),sqrt(1-eta.^2).*sin(az),eta];
[A,~]=far_fields(s,eta,az); radius=1e5*max(s.a,s.b); gap=0;
for j=[1,2*numel(s.theta)-1,2*numel(s.theta)]
    one=s; one.theta_inc=s.theta(ceil(j/2)); one.phi_inc=0;
    one.pol=double([mod(j,2)==1,mod(j,2)==0]);
    for k=1:numel(one.modes)
        for name={'a','b','c','d'}
            one.modes(k).(name{1})=s.modes(k).(name{1})(:,j);
        end
    end
    p=radius*dirs; p2=2*p;
    E1=spheroid_eval(one,p(:,1),p(:,2),p(:,3),'scattered');
    E2=spheroid_eval(one,p2(:,1),p2(:,2),p2(:,3),'scattered');
    extrapolated=2*(2*radius*exp(-2i*s.km*radius)*E2)-radius*exp(-1i*s.km*radius)*E1;
    gap=max(gap,relative(extrapolated,A(:,:,j)));
end
end

function run_sweep(s,out)
% Full baseline plus continuous and voxel convergence controls.
configs={ ...
    'baseline',256,.1,80,2,false; ...
    'window_51',512,.1,0,0,false; ...
    'sampling_0p05',1024,.05,0,0,false; ...
    'padding_3',256,.1,80,3,true; ...
    'voxel_0p0667',384,1/15,120,2,true; ...
    'voxel_window_51',512,.1,80,2,true};
allrows=table;
for k=1:size(configs,1)
    [name,N,dx,Nz,padding,subset]=configs{k,:};
    filename=fullfile(out,['oblate_',name,'.csv']);
    if subset
        ref=subset_reference(s,[1,6,11]); az=[0,45,90]*pi/180;
    else
        ref=s; az=(0:45:315)*pi/180;
    end
    rows=grid_sweep(ref,N,dx,Nz,padding,az,name);
    writetable(rows,filename);
    allrows=[allrows;rows]; %#ok<AGROW>
    writetable(allrows,fullfile(out,'oblate_angles.csv'));
end
fprintf('OBLATE_SWEEP_FINISHED rows=%d\n',height(allrows));
end

function s=subset_reference(s,indices)
s.theta=s.theta(indices); s.illumination_NA=s.illumination_NA(indices);
v=reshape([2*indices-1;2*indices],1,[]);
for k=1:numel(s.modes)
    for name={'a','b','c','d'}
        s.modes(k).(name{1})=s.modes(k).(name{1})(:,v);
    end
    s.modes(k).so=s.modes(k).so(:,indices); s.modes(k).si=s.modes(k).si(:,indices);
end
end

function results=grid_sweep(s,N,dx,Nz,padding,az,name)
timer=tic; k0=2*pi/s.lambda; k=s.km;
xy=(-N/2:N/2-1)*dx; [X,Y]=ndgrid(xy,xy);
freq=(-N/2:N/2-1)*2*pi/(N*dx); [KX,KY]=ndgrid(freq,freq);
prop=KX.^2+KY.^2<k^2; pupil=KX.^2+KY.^2<=(k0*s.NA_det)^2;
KZ=zeros(N); KZ(prop)=sqrt(k^2-KX(prop).^2-KY(prop).^2);
eta=KZ(pupil)/k; detector_phi=atan2(KY(pupil),KX(pupil));
factor=2*pi*1i*exp(1i*KZ(pupil)*s.z_det)./KZ(pupil);
% Each actual incident azimuth is evaluated on the same laboratory FFT grid.
thetas=[]; phis=[];
for phi=az
    ix=1:numel(s.theta); if phi~=0, ix=ix(s.theta(ix)>0); end
    thetas=[thetas;s.theta(ix)]; phis=[phis;repmat(phi,numel(ix),1)]; %#ok<AGROW>
end
use_voxel=Nz>0; volume_error=NaN;
if use_voxel
    z=(-Nz/2:Nz/2-1)*dx;
    inside=(xy(:).^2+xy.^2)/s.b^2+reshape(z.^2,1,1,[])/s.a^2<=1;
    n=s.n_m*ones(N,N,Nz); n(inside)=s.n_p;
    volume_error=nnz(inside)*dx^3/(4*pi*s.a*s.b^2/3)-1;
    fprintf('Voxel %s N=%d Nz=%d dx=%.6g padding=%d angles=%d\n', ...
        name,N,Nz,dx,padding,numel(thetas));
    [vB,vR,~,vi]=born_rytov_voxel(s.lambda,n,s.n_m,dx,thetas,phis, ...
        s.NA_det,s.z_det,struct('pad_factor',padding));
    clear n inside
    voxel_seconds=vi.time; clear vi
else
    voxel_seconds=NaN;
end
rows=struct([]); row=0;
for phi=az
    [EV,US]=far_fields(s,eta,detector_phi-phi);
    EV=EV.*factor; US=US.*factor;
    for j=1:numel(s.theta)
        if s.theta(j)==0 && phi~=0, continue; end
        theta=s.theta(j); ki=k*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
        eTM=[cos(theta),0,-sin(theta)]; eTE=[0,1,0];
        e0=cos(phi)*eTM-sin(phi)*eTE;
        % Compare in the incidence-plane frame; rotating all vectors back is unitary.
        exact=cos(phi)*EV(:,:,2*j-1)-sin(phi)*EV(:,:,2*j);
        co=exact*e0.'; scalar=US(:,j);
        [~,~,~,lab_e0]=incident_plane_wave(k,s.n_m,theta,phi,[1,0],0,0,0);
        rotation=[cos(phi),-sin(phi),0;sin(phi),cos(phi),0;0,0,1];
        assert(norm((rotation*e0.').'-lab_e0)<1e-13,'Polarization rotation mismatch.');
        [B,R,phase]=continuous_spectra(s,KX,KY,KZ,prop,X,Y,ki,dx);
        R(~pupil)=0;
        assert(all(isfinite(R),'all'),'Nonfinite continuous Rytov field.');
        row=row+1;
        rows(row).configuration=string(name);
        rows(row).Nxy=N; rows(row).Nz=Nz; rows(row).dx_um=dx;
        rows(row).padding=padding; rows(row).window_um=N*dx;
        rows(row).illumination_NA=s.illumination_NA(j);
        rows(row).theta_deg=theta*180/pi; rows(row).phi_deg=phi*180/pi;
        rows(row).continuous_born_scalar=relative(B(pupil),scalar);
        rows(row).continuous_rytov_scalar=relative(R(pupil),scalar);
        rows(row).continuous_born_co=relative(B(pupil),co);
        rows(row).continuous_rytov_co=relative(R(pupil),co);
        rows(row).continuous_born_full=relative(B(pupil)*e0,exact);
        rows(row).continuous_rytov_full=relative(R(pupil)*e0,exact);
        rows(row).scalar_co=relative(scalar,co);
        rows(row).scalar_full=relative(scalar*e0,exact);
        rows(row).polarization_leakage=relative(co*e0,exact);
        rows(row).voxel_born_scalar=NaN; rows(row).voxel_rytov_scalar=NaN;
        rows(row).voxel_born_co=NaN; rows(row).voxel_rytov_co=NaN;
        rows(row).voxel_born_full=NaN; rows(row).voxel_rytov_full=NaN;
        rows(row).voxel_born_continuous_gap=NaN; rows(row).voxel_rytov_continuous_gap=NaN;
        if use_voxel
            vb=ft(vB(:,:,row),dx); vr=ft(vR(:,:,row),dx);
            assert(norm(vr(~pupil))/norm(vr(:))<1e-12,'Rytov pupil leakage.');
            rows(row).voxel_born_scalar=relative(vb(pupil),scalar);
            rows(row).voxel_rytov_scalar=relative(vr(pupil),scalar);
            rows(row).voxel_born_co=relative(vb(pupil),co);
            rows(row).voxel_rytov_co=relative(vr(pupil),co);
            rows(row).voxel_born_full=relative(vb(pupil)*e0,exact);
            rows(row).voxel_rytov_full=relative(vr(pupil)*e0,exact);
            rows(row).voxel_born_continuous_gap=relative(vb(pupil),B(pupil));
            rows(row).voxel_rytov_continuous_gap=relative(vr(pupil),R(pupil));
        end
        dk=2*pi/(N*dx);
        rows(row).scalar_norm2=sum(abs(scalar).^2)*dk^2;
        rows(row).co_norm2=sum(abs(co).^2)*dk^2;
        rows(row).vector_norm2=sum(abs(exact).^2,'all')*dk^2;
        rows(row).max_abs_rytov_phase=max(abs(phase),[],'all');
        rows(row).voxel_volume_error=volume_error; rows(row).voxel_seconds=voxel_seconds;
        fprintf('%s NAi=%.2f phi=%3.0f continuous B/R co=%.3f/%.3f%% voxel=%.3f/%.3f%%\n', ...
            name,s.illumination_NA(j),phi*180/pi,100*rows(row).continuous_born_co, ...
            100*rows(row).continuous_rytov_co,100*rows(row).voxel_born_co, ...
            100*rows(row).voxel_rytov_co);
    end
end
results=struct2table(rows);
fprintf('GRID_FINISHED %s rows=%d elapsed=%.1fs\n',name,height(results),toc(timer));
end

function u=ift(S,dx)
u=fftshift(ifft2(ifftshift(S)))/dx^2;
end

function S=ft(u,dx)
S=fftshift(fft2(ifftshift(u)))*dx^2;
end

function summarize_results(out)
t=readtable(fullfile(out,'oblate_angles.csv'),'TextType','string');
names={'baseline','window_51','sampling_0p05'};
metrics={'continuous_born_scalar','continuous_rytov_scalar', ...
    'continuous_born_co','continuous_rytov_co','continuous_born_full', ...
    'continuous_rytov_full','voxel_born_scalar','voxel_rytov_scalar', ...
    'voxel_born_co','voxel_rytov_co','voxel_born_full','voxel_rytov_full', ...
    'scalar_co','scalar_full','polarization_leakage'};
rows=struct([]); n=0;
for k=1:numel(names)
    a=t(t.configuration==names{k},:); w=illumination_weights(a);
    assert(height(a)==81 && abs(sum(w)-1)<1e-14,'Angular coverage.');
    for m=1:numel(metrics)
        metric=metrics{m}; v=a.(metric); if any(isnan(v)), continue; end
        if endsWith(metric,'_scalar')
            norm2=a.scalar_norm2;
        elseif endsWith(metric,'_co')
            norm2=a.co_norm2;
        else
            norm2=a.vector_norm2;
        end
        n=n+1; rows(n).configuration=string(names{k}); rows(n).metric=string(metric);
        rows(n).pupil_weighted_mean=sum(w.*v);
        rows(n).pooled_L2=sqrt(sum(w.*v.^2.*norm2)/sum(w.*norm2));
        [rows(n).maximum,i]=max(v); rows(n).minimum=min(v);
        rows(n).worst_illumination_NA=a.illumination_NA(i);
        rows(n).worst_phi_deg=a.phi_deg(i);
    end
end
summary=struct2table(rows); writetable(summary,fullfile(out,'oblate_summary.csv'));
% Coarsening checks for the angular average; not a continuous-angle bound.
a=t(t.configuration=="sampling_0p05",:); w=illumination_weights(a);
fine_mean=sum(w.*a.continuous_rytov_co);
rad=a(abs(a.illumination_NA/.1-round(a.illumination_NA/.1))<1e-9,:);
az=a(abs(a.phi_deg/90-round(a.phi_deg/90))<1e-9,:);
radial_mean=sum(illumination_weights(rad).*rad.continuous_rytov_co);
azimuth_mean=sum(illumination_weights(az).*az.continuous_rytov_co);
angular_check=table(fine_mean,radial_mean,azimuth_mean, ...
    abs(radial_mean-fine_mean),abs(azimuth_mean-fine_mean),'VariableNames', ...
    {'fine_mean','radial_step_0p1_mean','azimuth_step_90_mean','radial_mean_gap','azimuth_mean_gap'});
writetable(angular_check,fullfile(out,'oblate_angular_convergence.csv'));
disp(summary(contains(summary.metric,'_co'),:));
disp(angular_check);
fprintf('OBLATE_SUMMARY_FINISHED\n');
end

function w=illumination_weights(t)
% Uniform area in the illumination pupil: annulus area, then equal azimuth weights.
q=unique(t.illumination_NA); edges=[0;(q(1:end-1)+q(2:end))/2;.5];
area=diff(edges.^2)/.5^2; w=zeros(height(t),1);
for j=1:numel(q)
    ix=t.illumination_NA==q(j); w(ix)=area(j)/nnz(ix);
end
end

function field_convergence(p,out)
% Compare spectra at identical physical frequencies, not just error magnitudes.
configs=[256,.1;512,.1;1024,.05;1024,.1]; rows=struct([]); row=0;
for ni=[0,.25,.5]
    az=[0,pi/4,pi/2]; if ni==0, az=0; end
    for phi=az
        theta=asin(ni/p.n_m); spectra=cell(4,1); pupils=cell(4,1);
        for j=1:4
            N=configs(j,1); dx=configs(j,2); k0=2*pi/p.lambda; k=k0*p.n_m;
            f=(-N/2:N/2-1)*2*pi/(N*dx); [KX,KY]=ndgrid(f,f);
            prop=KX.^2+KY.^2<k^2; pupil=KX.^2+KY.^2<=(k0*p.NA_det)^2;
            KZ=zeros(N); KZ(prop)=sqrt(k^2-KX(prop).^2-KY(prop).^2);
            ki=k*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
            x=(-N/2:N/2-1)*dx; [X,Y]=ndgrid(x,x);
            [~,R]=continuous_spectra(p,KX,KY,KZ,prop,X,Y,ki,dx);
            R(~pupil)=0;
            spectra{j}=R; pupils{j}=pupil;
        end
        row=row+1; rows(row).illumination_NA=ni; rows(row).phi_deg=phi*180/pi;
        a=spectra{1}; b=spectra{2}(1:2:end,1:2:end); mask=pupils{1};
        rows(row).window_25_to_51=relative(a(mask),b(mask));
        a=spectra{2}; b=spectra{3}(257:768,257:768); mask=pupils{2};
        rows(row).dx_0p1_to_0p05=relative(a(mask),b(mask));
        b=spectra{4}(1:2:end,1:2:end);
        rows(row).window_51_to_102=relative(a(mask),b(mask));
    end
end
results=struct2table(rows); writetable(results,fullfile(out,'oblate_field_convergence.csv'));
disp(results); assert(max(results.dx_0p1_to_0p05)<1e-5,'Continuous dx refinement.');
assert(max(results.window_51_to_102)<.005,'Continuous window refinement exceeds 0.5%%.');
fprintf('OBLATE_FIELD_CONVERGENCE_PASS\n');
end

function [B,R,phase]=continuous_spectra(p,KX,KY,KZ,prop,X,Y,ki,dx)
Q=sqrt(p.b^2*((KX(prop)-ki(1)).^2+(KY(prop)-ki(2)).^2) ...
    +p.a^2*(KZ(prop)-ki(3)).^2);
C=(2*pi/p.lambda)^2*(p.n_p^2-p.n_m^2)*4*pi*p.a*p.b^2/3;
B=complex(zeros(size(KX)));
B(prop)=1i*C*form_factor(Q).*exp(1i*KZ(prop)*p.z_det)./(2*KZ(prop));
u0=exp(1i*(ki(1)*X+ki(2)*Y+ki(3)*p.z_det));
phase=ift(B,dx)./u0; R=ft(u0.*expm1(phase),dx);
end

function validate_xz(s,out)
% Independent public vector evaluator and finite-difference scalar derivatives.
p=[0,0,0;.7,.2,.4;2,-.4,1;3.5,.3,0;0,0,5;1,.5,5;4,0,5;-6,.2,2];
[ex,u,du]=near_fields(s,p,true); public_gap=zeros(4,1);
for j=1:4
    one=s; one.theta_inc=s.theta(j); one.phi_inc=0; one.pol=[1,0];
    for k=1:numel(s.modes)
        for name={'a','b','c','d'}
            one.modes(k).(name{1})=s.modes(k).(name{1})(:,2*j-1);
        end
    end
    e=spheroid_eval(one,p(:,1),p(:,2),p(:,3),'total');
    public_gap(j)=relative(ex(:,j),e(:,1));
end
h=2e-5; fd=du; lap=zeros(size(u));
for dim=1:3
    pp=p; pm=p; pp(:,dim)=pp(:,dim)+h; pm(:,dim)=pm(:,dim)-h;
    [~,up,dp]=near_fields(s,pp,false); [~,um,dm]=near_fields(s,pm,false);
    fd(:,dim,:)=reshape((up-um)/(2*h),size(p,1),1,4);
    lap=lap+reshape((dp(:,dim,:)-dm(:,dim,:))/(2*h),size(u));
end
gradient_gap=relative(fd,du);
f=(2*pi/s.lambda)^2*(s.n_p^2-s.n_m^2)* ...
    ((p(:,1)/s.b).^2+(p(:,2)/s.b).^2+(p(:,3)/s.a).^2<1);
pde_gap=zeros(4,1);
for j=1:4
    ki=s.km*[sin(s.theta(j)),0,cos(s.theta(j))];
    loggrad=du(:,:,j)./u(:,j); dpsi=loggrad-1i*ki;
    lappsi=lap(:,j)./u(:,j)-sum(loggrad.^2,2);
    residual=lappsi+2i*sum(ki.*dpsi,2)+sum(dpsi.^2,2)+f;
    pde_gap(j)=norm(residual)/norm(f);
end
assert(max(public_gap)<1e-10 && gradient_gap<1e-6 && max(pde_gap)<1e-5, ...
    'Public vector and logarithmic PDE validation.');
t=table(max(public_gap),gradient_gap,max(pde_gap), ...
    'VariableNames',{'public_Ex_gap','gradient_fd_gap','log_pde_residual'});
writetable(t,fullfile(out,'oblate_xz_derivative_validation.csv'));
fprintf('XZ_PUBLIC_AND_PDE_PASS vector=%.3g gradient=%.3g PDE=%.3g\n', ...
    max(public_gap),gradient_gap,max(pde_gap));
end

function run_xz(s,out)
% Raw total Ex on y=0, not a detector-pupil-filtered section through the object.
x=linspace(-8,8,201); z=linspace(-5,9,176); [X,Z]=meshgrid(x,z);
points=[X(:),zeros(numel(X),1),Z(:)]; theta_deg=s.theta*180/pi;
u0=exp(1i*s.km*(points(:,1)*sin(s.theta.')+points(:,3)*cos(s.theta.')));
inside=(X(:)/s.b).^2+(Z(:)/s.a).^2<1-64*eps;
[exact,scalar,grad]=near_fields(s,points,true);
f0=(2*pi/s.lambda)^2*(s.n_p^2-s.n_m^2);
Q=complex(zeros(size(scalar)));
for j=1:numel(s.theta)
    dpsi=grad(:,:,j)./scalar(:,j)-1i*s.km*[sin(s.theta(j)),0,cos(s.theta(j))];
    Q(:,j)=sum(dpsi.^2,2);
end
save(fullfile(out,'oblate_xz_exact.mat'),'x','z','theta_deg','exact','scalar','grad','Q','u0','inside','f0','-v7');
% Validate raw near fields at independent points against a lower cutoff.
ids=unique(round(linspace(1,size(points,1),257)));
low=batch_reference(s.parameters,109,55);
[el,sl]=near_fields(low,points(ids,:),true);
near_vector_gap=relative(el,exact(ids,:)); near_scalar_gap=relative(sl,scalar(ids,:));
assert(max(near_vector_gap,near_scalar_gap)<1e-5,'Near-field cutoff check.');
clear low
% The derivative at zero scattering potential is the full-Green first Born field.
h=.005; plus=scalar_scaled_reference(s,h); minus=scalar_scaled_reference(s,-h);
born=complex(zeros(size(scalar)));
derivative=plus;
for k=1:numel(s.modes)
    derivative.modes(k).so=(plus.modes(k).so-minus.modes(k).so)/(2*h);
end
[~,bo]=near_fields(derivative,points(~inside,:),false);
born(~inside,:)=bo-u0(~inside,:);
[~,up]=near_fields(plus,points(inside,:),false);
[~,um]=near_fields(minus,points(inside,:),false);
born(inside,:)=(up-um)/(2*h);
rytov=u0.*exp(born./u0);
bornEx=(u0+born).*cos(s.theta.'); rytovEx=rytov.*cos(s.theta.');
save(fullfile(out,'oblate_xz_fields.mat'),'x','z','theta_deg','exact','bornEx','rytovEx', ...
    'scalar','u0','born','Q','inside','f0','near_vector_gap','near_scalar_gap','-v7');
fprintf('XZ_MAPS_SAVED near vector/scalar cutoff gaps %.3g %.3g\n',near_vector_gap,near_scalar_gap);

% Independent spatial Born check using the existing surface-integral implementation.
probe=[0,0,0;.7,0,.4;2,0,1;3.5,0,0;0,0,5;1,0,5;4,0,5;-6,0,2];
[~,bp]=near_fields(plus,probe,false); [~,bm]=near_fields(minus,probe,false);
born_probe=(bp-bm)/(2*h); born_surface_gap=zeros(numel(s.theta),1);
for j=1:numel(s.theta)
    [bs,~,~,info]=born_rytov_spheroid(s.lambda,s.n_p,s.n_m,s.a,s.b,s.theta(j),0, ...
        probe(:,1),probe(:,2),probe(:,3),struct('tol',1e-7));
    assert(info.validated); born_surface_gap(j)=relative(born_probe(:,j),bs);
end
assert(max(born_surface_gap)<1e-4,'Born derivative vs surface integral.');

% Detector diagnostics: extract the second logarithmic coefficient, without fitting it.
N=512; dx=.1; f=(-N/2:N/2-1)*2*pi/(N*dx); [KX,KY]=ndgrid(f,f);
prop=KX.^2+KY.^2<s.km^2; pupil=KX.^2+KY.^2<=(2*pi*s.NA_det/s.lambda)^2;
KZ=zeros(N); KZ(prop)=sqrt(s.km^2-KX(prop).^2-KY(prop).^2);
eta=KZ(prop)/s.km; phi=atan2(KY(prop),KX(prop));
factor=2*pi*1i*exp(1i*KZ(prop)*s.z_det)./KZ(prop);
[~,Uplus]=far_fields(plus,eta,phi); [~,Uminus]=far_fields(minus,eta,phi);
first=(Uplus-Uminus).*factor/(2*h); second=(Uplus+Uminus).*factor/(2*h^2);
clear plus minus derivative
plus2=scalar_scaled_reference(s,2*h); minus2=scalar_scaled_reference(s,-2*h);
[~,Up2]=far_fields(plus2,eta,phi); [~,Um2]=far_fields(minus2,eta,phi);
second_coarse=(Up2+Um2).*factor/(8*h^2);
second_step_gap=relative(second_coarse,second);
[~,bpp]=near_fields(plus2,probe,false); [~,bmm]=near_fields(minus2,probe,false);
born_step_gap=relative((bpp-bmm)/(4*h),born_probe);
assert(second_step_gap<1e-3 && born_step_gap<1e-4,'Potential-step refinement.');
clear plus2 minus2
[~,Uexact]=far_fields(s,eta,phi); Uexact=Uexact.*factor;
xy=(-N/2:N/2-1)*dx; [XX,YY]=ndgrid(xy,xy);
rows=struct([]); row=0; termrows=struct([]);
for j=1:numel(s.theta)
    ki=s.km*[sin(s.theta(j)),0,cos(s.theta(j))];
    [B,~,psi1]=continuous_spectra(s,KX,KY,KZ,prop,XX,YY,ki,dx);
    first_gap=relative(first(:,j),B(prop)); assert(first_gap<1e-4,'Weak derivative vs analytic Born.');
    inc=exp(1i*(ki(1)*XX+ki(3)*s.z_det));
    U2=zeros(N); U2(prop)=second(:,j);
    psi2=ift(U2,dx)./inc-.5*psi1.^2;
    for t=[.1,.25,.5,1]
        if t==1
            target=Uexact(:,j);
        else
            filename=fullfile(out,sprintf('oblate_contrast_t%g.mat',t));
            if j==1
                scaled=scalar_scaled_reference(s,t); [~,Ut]=far_fields(scaled,eta,phi);
                Ut=Ut.*factor; save(filename,'Ut','-v7'); clear scaled
            else
                loaded=load(filename,'Ut'); Ut=loaded.Ut;
            end
            target=Ut(:,j);
        end
        R1=ft(inc.*expm1(t*psi1),dx); R2=ft(inc.*expm1(t*psi1+t^2*psi2),dx);
        target_grid=zeros(N); target_grid(prop)=target;
        row=row+1; rows(row).theta_deg=theta_deg(j); rows(row).potential_scale=t;
        rows(row).index=sqrt(s.n_m^2+t*(s.n_p^2-s.n_m^2));
        rows(row).born_error=relative(t*B(pupil),target_grid(pupil));
        rows(row).rytov1_error=relative(R1(pupil),target_grid(pupil));
        rows(row).rytov2_error=relative(R2(pupil),target_grid(pupil));
        rows(row).psi2_over_psi1=norm(t^2*psi2(:))/norm(t*psi1(:));
        fprintf('OBLATE_CAUSE theta=%g t=%.2f R1=%.4f%% R2=%.4f%%\n',theta_deg(j),t, ...
            100*rows(row).rytov1_error,100*rows(row).rytov2_error);
    end
    valid=inside & abs(scalar(:,j))>.1;
    termrows(j).theta_deg=theta_deg(j);
    termrows(j).slice_Q_over_f_rms=norm(Q(valid,j))/sqrt(nnz(valid))/f0;
    termrows(j).minimum_interior_abs_U=min(abs(scalar(inside,j)));
    termrows(j).masked_interior_fraction=1-nnz(valid)/nnz(inside);
    termrows(j).xz_total_Ex_born_error=relative(bornEx(:,j),exact(:,j));
    termrows(j).xz_total_Ex_rytov_error=relative(rytovEx(:,j),exact(:,j));
    termrows(j).xz_total_Ex_scalar_error=relative(scalar(:,j)*cos(s.theta(j)),exact(:,j));
end
writetable(struct2table(rows),fullfile(out,'oblate_rytov_cause.csv'));
writetable(struct2table(termrows),fullfile(out,'oblate_xz_metrics.csv'));
validation=table(near_vector_gap,near_scalar_gap,max(born_surface_gap),born_step_gap,second_step_gap, ...
    'VariableNames',{'near_vector_cutoff','near_scalar_cutoff','born_surface_gap','born_step_gap','second_step_gap'});
writetable(validation,fullfile(out,'oblate_xz_validation.csv'));
fprintf('OBLATE_XZ_AND_CAUSE_PASS\n');
end

function s=scalar_scaled_reference(base,t)
% Same geometry/exterior basis; only scalar boundary coefficients are replaced.
s=base; s.scalar_only=true;
s.n_p=sqrt(s.n_m^2+t*(base.parameters.n_p^2-s.n_m^2));
s.c_int=2*pi*s.n_p/s.lambda*s.semifocal;
[g,w]=meridian(s,s.Qtheta); [gv,wv]=meridian(s,s.Qtheta+17);
inc=boundary_incident(s,g); iv=boundary_incident(s,gv);
n=numel(s.theta); err=zeros(1,n); timer=tic;
for mu=0:s.M
    d=s.mu_data{mu+1}; d.eg_int=sph_eigen('oblate',mu,s.c_int,s.L);
    s.mu_data{mu+1}=d;
    Wo=sph_angular(d.eg_ext,g); Wi=sph_angular(d.eg_int,g);
    Vo=sph_radial_ode_eval(d.rod_ext,g); [Vi,ok]=sph_radial(d.eg_int,1,g,s.c_int); assert(ok);
    Wov=sph_angular(d.eg_ext,gv); Wiv=sph_angular(d.eg_int,gv);
    Vov=sph_radial_ode_eval(d.rod_ext,gv); [Viv,ok]=sph_radial(d.eg_int,1,gv,s.c_int); assert(ok);
    A=scalar_matrix(Wo,Vo,Wi,Vi,g,s.c_ext); Av=scalar_matrix(Wov,Vov,Wiv,Viv,gv,s.c_ext);
    for k=find([s.modes.mu]==mu)
        bin=mod(s.modes(k).m,s.Qphi)+1;
        rhs=-reshape(inc.s(:,bin,:),[],n);
        c=scaled_solve(A.*sqrt(repmat(w,2,1)),rhs.*sqrt(repmat(w,2,1)));
        count=s.L-mu+1; s.modes(k).so=c(1:count,:); s.modes(k).si=c(count+1:end,:);
        residual=Av*c+reshape(iv.s(:,bin,:),[],n);
        err=err+sum(repmat(wv,2,1).*abs(residual).^2,1);
    end
end
den=squeeze(sum(sum(repmat(wv,2,1).*abs(iv.s).^2,1),2)).';
assert(max(sqrt(err./den))<1e-7,'Scaled scalar boundary residual.');
fprintf('SCALED_REFERENCE t=%g residual=%.3g elapsed=%.1fs\n',t,max(sqrt(err./den)),toc(timer));
end

function [Ex,U,dU]=near_fields(s,points,vector)
% Batched total field and optional exact Cartesian scalar gradient.
n=numel(s.theta); count=size(points,1); gradient=nargout>2;
if vector, assert(~isfield(s,'scalar_only') || ~s.scalar_only); end
Ex=complex(zeros(count,n)); U=Ex;
if gradient, dU=complex(zeros(count,3,n)); else, dU=[]; end
for first=1:1000:count
    ix=first:min(first+999,count); p=points(ix,:);
    outer=(p(:,1)/s.b).^2+(p(:,2)/s.b).^2+(p(:,3)/s.a).^2>=1-64*eps;
    for exterior=[false,true]
        jj=find(outer==exterior); if isempty(jj), continue; end
        ids=ix(jj); q=p(jj,:);
        g=sph_coords('oblate',q(:,1)/s.semifocal,q(:,2)/s.semifocal,q(:,3)/s.semifocal);
        for mu=0:s.M
            data=s.mu_data{mu+1};
            if exterior
                W=sph_angular(data.eg_ext,g); V=sph_radial_ode_eval(data.rod_ext,g); c=s.c_ext;
            else
                W=sph_angular(data.eg_int,g); [V,ok]=sph_radial(data.eg_int,1,g,s.c_int); assert(ok); c=s.c_int;
            end
            S=(W.W0.*V.V0).';
            if gradient
                De=(W.W1.*V.V0).'./(s.semifocal*sqrt(g.D));
                Dxi=(W.W0.*V.V1).'./(s.semifocal*sqrt(g.D));
                Dphi=(W.Wm.*(V.V0./sqrt(g.alpha.'))).'/s.semifocal;
            end
            for k=find([s.modes.mu]==mu)
                mode=s.modes(k); phase=exp(1i*mode.m*g.phi);
                if exterior, coeff=mode.so; else, coeff=mode.si; end
                U(ids,:)=U(ids,:)+phase.*(S*coeff);
                if gradient
                    de=phase.*(De*coeff); dxi=phase.*(Dxi*coeff);
                    dp=1i*sign(mode.m)*phase.*(Dphi*coeff);
                    for dim=1:3
                        dU(ids,dim,:)=dU(ids,dim,:)+reshape(de.*g.e_eta(:,dim)+ ...
                            dxi.*g.e_xi(:,dim)+dp.*g.e_phi(:,dim),numel(ids),1,n);
                    end
                end
                if vector
                    [Mw,Nw]=sph_vecwave(W,V,g,mode.m,c);
                    if exterior, ca=mode.a(:,1:2:end); cb=mode.b(:,1:2:end);
                    else, ca=mode.c(:,1:2:end); cb=mode.d(:,1:2:end); end
                    Ex(ids,:)=Ex(ids,:)+phase.*((Nw.eta.'*ca+Mw.eta.'*cb).*g.e_eta(:,1)+ ...
                        (Nw.xi.'*ca+Mw.xi.'*cb).*g.e_xi(:,1)+ ...
                        (Nw.phi.'*ca+Mw.phi.'*cb).*g.e_phi(:,1));
                end
            end
        end
        if exterior
            inc=exp(1i*s.km*(q(:,1)*sin(s.theta.')+q(:,3)*cos(s.theta.')));
            U(ids,:)=U(ids,:)+inc;
            if vector, Ex(ids,:)=Ex(ids,:)+inc.*cos(s.theta.'); end
            if gradient
                dU(ids,1,:)=dU(ids,1,:)+reshape(1i*s.km*inc.*sin(s.theta.'),numel(ids),1,n);
                dU(ids,3,:)=dU(ids,3,:)+reshape(1i*s.km*inc.*cos(s.theta.'),numel(ids),1,n);
            end
        end
    end
    if mod(first-1,5000)==0, fprintf('NEAR_FIELDS %d/%d vector=%d gradient=%d\n',min(first+999,count),count,vector,gradient); end
end
assert(all(isfinite(U),'all') && all(isfinite(Ex),'all'),'Nonfinite near field.');
end

function f=form_factor(Q)
f=1-Q.^2/10+Q.^4/280; ii=abs(Q)>1e-3; q=Q(ii);
f(ii)=3*(sin(q)-q.*cos(q))./q.^3;
end

function e=relative(a,b)
e=norm(a(:)-b(:))/norm(b(:));
end
