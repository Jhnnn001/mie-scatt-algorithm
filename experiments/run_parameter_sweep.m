function run_parameter_sweep(project,out)
%RUN_PARAMETER_SWEEP Shared, checkpointed size/contrast study. All lengths in um.
assert(ismember(project,{'size','contrast'}),'Unknown experiment.');
root=fileparts(fileparts(out));
addpath(fullfile(root,'experiments','field-analysis'),fullfile(root,'maxwell-solver'));
api=oblate_field_analysis('helpers');
oblate_field_analysis('check');
p=struct('lambda',.532,'n_m',1.335381534,'n_p',1.365,'a',3,'b',5, ...
    'NA_det',1,'z_det',5,'illumination_NA',(0:.05:.5).');
p.theta=asin(p.illumination_NA/p.n_m);
if strcmp(project,'size'), factors=[1,.1,.01]; else, factors=[1,.5,.01]; end
allrows=table; allchecks=table; configs=struct([]);
for ic=1:3
    q=p; factor=factors(ic); name=sprintf('%s_%g_percent',project,100*factor);
    if strcmp(project,'size')
        q.a=p.a*factor; q.b=p.b*factor;
        orders=[114,59;36,18;40,14]; L=orders(ic,1); M=orders(ic,2);
    else
        q.n_p=p.n_m+factor*(p.n_p-p.n_m); L=114; M=59;
    end
    config=struct('case_id',ic,'factor',factor,'a_um',q.a,'b_um',q.b, ...
        'n_particle',q.n_p,'n_background',q.n_m,'index_ratio',q.n_p/q.n_m, ...
        'delta_n',q.n_p-q.n_m,'lambda0_um',q.lambda,'lambda_medium_um',q.lambda/q.n_m, ...
        'phase_diameter_rad',2*pi/q.lambda*(q.n_p-q.n_m)*2*q.a, ...
        'kb',2*pi*q.n_m/q.lambda*q.b,'volume_um3',4*pi*q.a*q.b^2/3, ...
        'z_detector_um',q.z_det,'L',L,'M',M);
    if ic==1, configs=config; else, configs(ic)=config; end
    writetable(struct2table(configs),fullfile(out,'parameters.csv'));
    checkpoint=fullfile(out,[name,'_reference.mat']);
    validation=fullfile(out,[name,'_reference_validation.csv']);
    if isfile(checkpoint) && isfile(validation)
        load(checkpoint,'s');
        assert(isequal(s.parameters,q) && s.L==L && s.M==M,'Checkpoint mismatch.');
        check=readtable(validation);
    else
        if factor==1
            high=load(fullfile(root,'experiments','field-analysis','oblate_reference_fine.mat'),'s'); s=high.s;
            low=load(fullfile(root,'experiments','field-analysis','oblate_reference_coarse.mat'),'s'); low=low.s;
            assert(isequal(s.parameters,q),'Original reference parameters differ.');
            public_gap=s.public_far_gap;
        else
            s=api.reference(q,L,M);
            low=api.reference(q,L-5,M-4);
            public_gap=api.public_far_check(s);
        end
        [g,w]=gauss_legendre(48); na=(g+1)/2; w=w/2;
        [NA,PHI]=ndgrid(na,2*pi*(0:63)/64);
        eta=sqrt(1-(NA/q.n_m).^2);
        [e,u]=api.far_fields(s,eta(:),PHI(:));
        [el,ul]=api.far_fields(low,eta(:),PHI(:));
        weights=repmat(w.*na./(q.n_m^2-na.^2),64,1);
        vg=rel(el.*sqrt(weights),e.*sqrt(weights));
        sg=rel(ul.*sqrt(weights),u.*sqrt(weights));
        assert(max([vg,sg])<1e-6 && public_gap<1e-5,'Reference convergence.');
        check=table(ic,vg,sg,s.vector_boundary_error,s.scalar_boundary_error,public_gap, ...
            'VariableNames',{'case_id','cutoff_vector','cutoff_scalar', ...
            'boundary_vector','boundary_scalar','public_far_gap'});
        % Detector spectra need coefficients, not dense outgoing trajectories.
        for mu=1:numel(s.mu_data)
            if isfield(s.mu_data{mu}.rod_ext,'solution')
                s.mu_data{mu}.rod_ext=rmfield(s.mu_data{mu}.rod_ext,'solution');
            end
        end
        save(checkpoint,'s','-v7'); writetable(check,validation);
        clear low high el ul e u
    end
    allchecks=[allchecks;check]; %#ok<AGROW>
    writetable(allchecks,fullfile(out,'reference_validation.csv'));
    output=fullfile(out,[name,'_angles.csv']);
    if isfile(output)
        rows=readtable(output,'TextType','string');
        assert(height(rows)==81 && all(rows.case_id==ic),'Incomplete angle checkpoint.');
    else
        rows=detector_sweep(api,s,512,.1,(0:45:315)*pi/180,ic);
        writetable(rows,output);
    end
    allrows=[allrows;rows]; %#ok<AGROW>
    writetable(allrows,fullfile(out,'angles.csv'));
    convergence_file=fullfile(out,[name,'_detector_convergence.csv']);
    if ~isfile(convergence_file)
        small=api.subset(s,[1,6,11]); az=[0,pi/4,pi/2];
        [coarse,fields]=detector_sweep(api,small,512,.1,az,ic);
        [fine,fields_fine]=detector_sweep(api,small,1024,.05,az,ic);
        [wide,fields_wide]=detector_sweep(api,small,1024,.1,az,ic);
        sample_gap=zeros(height(coarse),1); window_gap=sample_gap;
        for j=1:height(coarse)
            a=fields{j}; b=fields_fine{j}(257:768,257:768);
            c=fields_wide{j}(1:2:end,1:2:end);
            sample_gap(j)=rel(a,b); window_gap(j)=rel(a,c);
        end
        convergence=table(repmat(ic,height(coarse),1),coarse.illumination_NA,coarse.phi_deg, ...
            sample_gap,window_gap, ...
            abs(coarse.rytov_scalar-fine.rytov_scalar), ...
            abs(coarse.rytov_scalar-wide.rytov_scalar), ...
            abs(coarse.rytov_co-wide.rytov_co), ...
            'VariableNames',{'case_id','illumination_NA','phi_deg','sampling_field_gap', ...
            'window_field_gap','sampling_error_gap','window_error_gap','window_co_error_gap'});
        writetable(convergence,convergence_file);
        clear fields fields_fine fields_wide
    end
    convergence=readtable(convergence_file);
    assert(max(convergence.sampling_field_gap)<1e-4,'Detector sampling refinement.');
    assert(max(convergence.window_field_gap)<.002,'Detector window refinement exceeds 0.2 percent.');
    fprintf('CASE_COMPLETE project=%s case=%d factor=%.8g\n',project,ic,factor);
end
summary=summarize(api,allrows); writetable(summary,fullfile(out,'summary.csv'));
assert(height(allrows)==243 && all(isfinite(summary.mean_error)),'Complete finite dataset.');
disp(summary(ismember(summary.metric,["born_scalar","rytov_scalar","born_co","rytov_co","vector_born_full"]),:));
fprintf('PARAMETER_EXPERIMENT_PASS project=%s angles=%d\n',project,height(allrows));
end

function [rows,fields]=detector_sweep(api,s,N,dx,az,ic)
k=s.km; k0=2*pi/s.lambda;
x=(-N/2:N/2-1)*dx; [X,Y]=ndgrid(x,x);
f=(-N/2:N/2-1)*2*pi/(N*dx); [KX,KY]=ndgrid(f,f);
prop=KX.^2+KY.^2<k^2; pupil=KX.^2+KY.^2<=(k0*s.NA_det)^2;
KZ=zeros(N); KZ(prop)=sqrt(k^2-KX(prop).^2-KY(prop).^2);
eta=KZ(pupil)/k; detector_phi=atan2(KY(pupil),KX(pupil));
factor=2*pi*1i*exp(1i*KZ(pupil)*s.z_det)./KZ(pupil);
row=0; data=struct([]); fields={};
for phi=az
    [EV,US]=api.far_fields(s,eta,detector_phi-phi);
    EV=EV.*factor; US=US.*factor;
    for j=1:numel(s.theta)
        theta=s.theta(j); if theta==0 && phi~=0, continue; end
        ki=k*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
        e0=cos(phi)*[cos(theta),0,-sin(theta)]-sin(phi)*[0,1,0];
        exact=cos(phi)*EV(:,:,2*j-1)-sin(phi)*EV(:,:,2*j);
        co=exact*e0.'; scalar=US(:,j);
        [~,~,~,lab_e0]=incident_plane_wave(k,s.n_m,theta,phi,[1,0],0,0,0);
        rot=[cos(phi),-sin(phi),0;sin(phi),cos(phi),0;0,0,1];
        assert(norm((rot*e0.').'-lab_e0)<1e-13,'Polarization mismatch.');
        [B,R,phase]=api.spectra(s,KX,KY,KZ,prop,X,Y,ki,dx);
        R(~pupil)=0; B=B(pupil); r=R(pupil);
        % Vector Born is the same Born amplitude with the Maxwell projector.
        detdir=[sqrt(1-eta.^2).*cos(detector_phi-phi), ...
            sqrt(1-eta.^2).*sin(detector_phi-phi),eta];
        VB=B.*(e0-(detdir*e0.').*detdir);
        row=row+1;
        data(row).case_id=ic; data(row).illumination_NA=s.illumination_NA(j);
        data(row).theta_deg=theta*180/pi; data(row).phi_deg=phi*180/pi;
        data(row).born_scalar=rel(B,scalar); data(row).rytov_scalar=rel(r,scalar);
        data(row).born_co=rel(B,co); data(row).rytov_co=rel(r,co);
        data(row).born_full=rel(B*e0,exact); data(row).rytov_full=rel(r*e0,exact);
        data(row).scalar_co=rel(scalar,co); data(row).scalar_full=rel(scalar*e0,exact);
        data(row).vector_born_co=rel(VB*e0.',co);
        data(row).vector_born_full=rel(VB,exact);
        data(row).scalar_norm2=sum(abs(scalar).^2)*(2*pi/(N*dx))^2;
        data(row).co_norm2=sum(abs(co).^2)*(2*pi/(N*dx))^2;
        data(row).vector_norm2=sum(abs(exact).^2,'all')*(2*pi/(N*dx))^2;
        data(row).max_abs_rytov_phase=max(abs(phase),[],'all');
        if nargout>1, fields{row,1}=R; end %#ok<AGROW>
        fprintf('DETECTOR case=%d N=%d dx=%.3g NAi=%.2f phi=%g scalar B/R=%.5f/%.5f%% co=%.5f/%.5f%%\n', ...
            ic,N,dx,s.illumination_NA(j),phi*180/pi,100*data(row).born_scalar, ...
            100*data(row).rytov_scalar,100*data(row).born_co,100*data(row).rytov_co);
    end
end
rows=struct2table(data);
assert(all(isfinite(rows{:,5:end}),'all'),'Nonfinite detector results.');
end

function result=summarize(api,rows)
metrics={'born_scalar','rytov_scalar','born_co','rytov_co','born_full','rytov_full', ...
    'scalar_co','scalar_full','vector_born_co','vector_born_full'};
data=struct([]); k=0;
for ic=1:3
    t=rows(rows.case_id==ic,:); w=api.weights(t);
    assert(height(t)==81 && abs(sum(w)-1)<1e-14,'Illumination coverage.');
    for j=1:numel(metrics)
        metric=metrics{j}; v=t.(metric);
        if endsWith(metric,'scalar'), norm2=t.scalar_norm2;
        elseif endsWith(metric,'co'), norm2=t.co_norm2; else, norm2=t.vector_norm2; end
        k=k+1; data(k).case_id=ic; data(k).metric=string(metric);
        data(k).mean_error=sum(w.*v);
        data(k).pooled_L2=sqrt(sum(w.*v.^2.*norm2)/sum(w.*norm2));
        data(k).minimum=min(v); [data(k).maximum,ix]=max(v);
        data(k).worst_NA=t.illumination_NA(ix); data(k).worst_phi=t.phi_deg(ix);
        data(k).fraction_below_1pct=sum(w(v<.01)); data(k).fraction_below_5pct=sum(w(v<.05));
        coarse=t(abs(t.illumination_NA/.1-round(t.illumination_NA/.1))<1e-9,:);
        data(k).radial_coarsening_gap=abs(sum(api.weights(coarse).*coarse.(metric))-sum(w.*v));
        coarse=t(abs(t.phi_deg/90-round(t.phi_deg/90))<1e-9,:);
        data(k).azimuth_coarsening_gap=abs(sum(api.weights(coarse).*coarse.(metric))-sum(w.*v));
    end
end
result=struct2table(data);
end

function e=rel(a,b)
den=norm(b(:)); assert(den>0 && isfinite(den),'Reference norm is not positive.');
e=norm(a(:)-b(:))/den;
end
