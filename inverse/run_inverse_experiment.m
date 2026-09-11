function summary = run_inverse_experiment(name, opts)
%RUN_INVERSE_EXPERIMENT Saved, independently generated low-contrast ODT study.
% Examples: run_inverse_experiment('baseline');
% run_inverse_experiment('matched',struct('source','matched'));
% The true support is used ONLY for error metrics and the labeled control.
if nargin<1, name='baseline'; end
if nargin<2, opts=struct; end
assert(ischar(name) && ~isempty(regexp(name,'^[a-z0-9_]+$','once')),'Use a simple run label.');
defaults=struct('source','exact','grid_size',[128,128,80],'dx',.1,'padding',2, ...
    'gp_iter',40,'tv_inner',100,'outer_iter',5,'alpha_relative',4, ...
    'tol',1e-4,'cg_tol',1e-4,'cg_max_iter',50,'discrepancy',.01, ...
    'methods',{{'direct','gp','tv'}},'implementation',3,'rho',.1);
assert(isstruct(opts) && isscalar(opts) && ...
    all(ismember(fieldnames(opts),fieldnames(defaults))),'Unknown experiment option.');
cfg=defaults;
for key=fieldnames(opts).', cfg.(key{1})=opts.(key{1}); end
cfg.source=validatestring(cfg.source,{'exact','matched'});
validateattributes(cfg.rho,{'numeric'},{'real','finite','scalar','positive'});
assert(all(ismember(cfg.methods,{'direct','gp','tv'})),'Unknown reconstruction method.');
here=fileparts(mfilename('fullpath')); addpath(here);
out=fullfile(here,'results'); if ~isfolder(out), mkdir(out); end
% Validate a reused label before touching any of its report artifacts.
for prior={'geometry','direct','gp','tv'}
    filename=fullfile(out,[name,'_',prior{1},'.mat']);
    if isfile(filename)
        if strcmp(prior{1},'geometry'), saved=load(filename,'cfg');
        else, saved=load(filename,'R'); end
        previous=cfg;
        if isfield(saved,'R'), previous=saved.R.config; end
        if isfield(saved,'cfg'), previous=saved.cfg; end
        if ~isfield(previous,'rho'), previous.rho=.1; end
        assert(isequal(orderfields(previous),orderfields(cfg)),'Existing run settings differ; use a new label.');
    end
end
logfile=fullfile(out,[name,'.log']); diary(logfile); cleanup=onCleanup(@() diary('off'));
fprintf('INVERSE_RUN %s source=%s padding=%g grid=%s dx=%g\n', ...
    name,cfg.source,cfg.padding,mat2str(cfg.grid_size),cfg.dx);
cache=fullfile(out,'exact_detector_512.mat');
rebuild_cache=~isfile(cache);
if ~rebuild_cache
    loaded=load(cache,'raw'); raw=loaded.raw;
    p=raw.parameters; gg=raw.geometry;
    assert(abs(p.n_p-1.33567771866)<1e-12 && p.n_m==1.335381534 && ...
        p.lambda==.532 && p.a==3 && p.b==5 && p.z_det==5 && p.NA_det==1 && ...
        gg.NA_det==1 && gg.z_det==5 && gg.detector_dx==.1 && ...
        numel(gg.kx)==512 && numel(gg.ky)==512 && size(gg.k_incident,1)==81 && ...
        size(raw.co_spectrum,2)==81 && numel(raw.samples)==81*nnz(gg.pupil), ...
        'Wrong measurement cache.');
    expected_angles=readtable(fullfile(here,'..','spheroid-analytic-forward','2','case_3_angles.csv'));
    assert(isequal(raw.angles,expected_angles),'Cached illumination/order differs from the validated study.');
    rebuild_cache=~isfield(raw,'preparation_version');
    if ~rebuild_cache, assert(raw.preparation_version==2,'Unsupported measurement preparation version.'); end
end
if rebuild_cache
    previous_samples=[];
    if exist('raw','var'), previous_samples=raw.samples; end
    m=odt_exact_data;
    data=odt_prepare(m.U,m.u0,m.geometry,'rytov',struct('weights',m.weights));
    raw=struct('q',data.q,'samples',data.samples,'geometry',m.geometry, ...
        'weights',m.weights,'angles',m.angles,'parameters',m.parameters, ...
        'co_spectrum',m.co_spectrum,'scalar_spectrum',m.scalar_spectrum, ...
        'phase_order_error',data.phase_order_error,'border_phase',data.border_phase, ...
        'reference_validation',m.reference_validation,'preparation_version',2);
    if ~isempty(previous_samples)
        assert(norm(previous_samples-raw.samples)/norm(raw.samples)<1e-12, ...
            'Unversioned measurement cache did not reproduce; preserve it and use a separate result folder.');
    end
    save(cache,'raw','-v7'); clear m data
end
g=raw.geometry; g.grid_size=cfg.grid_size; g.dx=cfg.dx; g.pad_factor=cfg.padding;
assert((g.grid_size(3)/2-.5)*g.dx<g.z_det,'Detector must be above reconstruction volume.');
ns=nnz(g.pupil); weights=repelem(raw.weights,ns);
op=odt_operator(raw.q,cfg.grid_size,cfg.dx,cfg.padding,weights);
physical_y=raw.samples.*op.sample_scale;
nm=raw.parameters.n_m; np=raw.parameters.n_p; sz=cfg.grid_size;
x=(-sz(1)/2:sz(1)/2-1)*cfg.dx; y=(-sz(2)/2:sz(2)/2-1)*cfg.dx;
z=(-sz(3)/2:sz(3)/2-1)*cfg.dx; [X,Y,Z]=ndgrid(x,y,z);
inside=(X/5).^2+(Y/5).^2+(Z/3).^2<1;
truth=nm*ones(sz); truth(inside)=np; true_chi=truth.^2/nm^2-1;
clear X Y Z
truth_prediction=op.forward(true_chi);
truth_mismatch=norm(truth_prediction-physical_y)/norm(physical_y);
fprintf('Known object vs physical Rytov Ewald data: %.6f%%\n',100*truth_mismatch);
target=physical_y;
if strcmp(cfg.source,'matched'), target=truth_prediction; end
data=struct('operator',op,'y',target,'n_m',nm);
initial=op.direct(target);
% Keep the physical TV/data balance fixed when padding or voxel pitch changes.
% A scales physical Fourier data by c; physical TV is dx^2*sum(|D chi|).
alpha=cfg.alpha_relative*max(abs(initial(:)))*op.physical_normalization*cfg.dx^2;
physical_normalization=op.physical_normalization;
fprintf('Data-derived alpha=%.9g; %d Ewald samples.\n',alpha,numel(target));
% The mask visualizes interpolation support, including inferred Hermitian data.
[~,coverage]=op.gridding(target);
coverage_xz=squeeze(coverage(:,size(coverage,2)/2+1,:));
coverage_xy=coverage(:,:,size(coverage,3)/2+1); clear coverage
axis_qx=(-sz(1)/2:sz(1)/2-1)*2*pi/(sz(1)*cfg.dx);
axis_qz=(-sz(3)/2:sz(3)/2-1)*2*pi/(sz(3)*cfg.dx);
geometry_file=fullfile(out,[name,'_geometry.mat']);
save(geometry_file,'x','y','z','truth','coverage_xz','coverage_xy','axis_qx','axis_qz', ...
    'cfg','truth_mismatch','alpha','physical_normalization','nm','np','-v7');
rows=struct([]);
for method_cell=cfg.methods
    method=method_cell{1}; filename=fullfile(out,[name,'_',method,'.mat']);
    if isfile(filename)
        loaded=load(filename,'R'); R=loaded.R;
        if ~isfield(R.config,'rho'), R.config.rho=.1; end
        assert(isequal(orderfields(R.config),orderfields(cfg)),'Existing result settings differ; use a new run label.');
        fprintf('Reusing validated settings checkpoint %s\n',filename);
    else
        iterations=cfg.tv_inner;
        if strcmp(method,'gp'), iterations=cfg.gp_iter; end
        solver=struct('alpha',alpha,'max_iter',iterations,'outer_iter',cfg.outer_iter, ...
            'rho',cfg.rho,'rho_positive',cfg.rho, ...
            'tol',cfg.tol,'cg_tol',cfg.cg_tol,'cg_max_iter',cfg.cg_max_iter, ...
            'discrepancy',cfg.discrepancy,'verbose',true);
        [n,chi,info]=odt_reconstruct(data,method,solver);
        metric=metrics(n,truth,inside,nm,np,cfg.dx);
        metric.method=string(method); metric.source=string(cfg.source);
        metric.ewald_residual=info.relative_data_residual;
        metric.physical_ewald_residual=norm(op.forward(chi)-physical_y)/norm(physical_y);
        [~,~,prediction]=odt_predict(chi,g,cfg.padding);
        reference=raw.co_spectrum.*sqrt(raw.weights.');
        predicted=prediction.rytov.*sqrt(raw.weights.');
        metric.physical_scattered_field_error=norm(predicted(:)-reference(:))/norm(reference(:));
        metric.converged=info.converged; metric.inner_converged=info.inner_converged;
        metric.iterations=info.iterations; metric.outer_iterations=info.outer_iterations;
        metric.stationarity=NaN;
        if strcmp(method,'tv'), metric.stationarity=info.stationarity_residual; end
        metric.seconds=info.time;
        R=struct('n',n,'chi',chi,'info',info,'metric',metric,'config',cfg);
        save(filename,'R','-v7');
        history=array2table(info.history,'VariableNames', ...
            {'outer','inner','objective','data_residual','primal_or_step','dual','cg_residual','cg_iterations'});
        writetable(history,fullfile(out,[name,'_',method,'_history.csv']));
    end
    if ~strcmp(method,'tv') && ~isnan(R.metric.stationarity)
        R.metric.stationarity=NaN; R.info.stationarity_residual=NaN;
        save(filename,'R','-v7');
    end
    % Refresh evaluation-only metadata when reusing a numerical checkpoint.
    evaluated=metrics(R.n,truth,inside,nm,np,cfg.dx);
    for key=fieldnames(evaluated).', R.metric.(key{1})=evaluated.(key{1}); end
    R.metric=orderfields(R.metric);
    save(filename,'R','-v7');
    rows=[rows;R.metric]; %#ok<AGROW>
    summary=struct2table(rows);
    writetable(summary,fullfile(out,[name,'_metrics.csv']));
    fprintf('%s COMPLETE RI=%.5f%% field=%.5f%% data=%.5f%% converged=%d\n', ...
        method,100*R.metric.relative_ri_error,100*R.metric.physical_scattered_field_error, ...
        100*R.metric.ewald_residual,R.metric.converged);
    clear R n chi prediction predicted reference
end
assert(height(summary)==numel(cfg.methods) && all(isfinite(summary.relative_ri_error)), ...
    'Incomplete or nonfinite reconstruction results.');
disp(summary);
fprintf('INVERSE_EXPERIMENT_PASS %s\n',name);
end

function m=metrics(n,truth,inside,nm,np,dx)
contrast=np-nm; denominator=norm(truth(:)-nm);
m=struct;
m.relative_ri_error=norm(n(:)-truth(:))/denominator;
m.object_error=norm(n(inside)-np)/denominator;
m.background_error=norm(n(~inside)-nm)/denominator;
m.mean_object_n=mean(n(inside));
m.relative_index_bias=(m.mean_object_n-np)/contrast;
sz=size(n); profile=squeeze(n(sz(1)/2+1,sz(2)/2+1,:));
truth_profile=squeeze(truth(sz(1)/2+1,sz(2)/2+1,:));
above=find(profile>nm+.5*contrast); true_above=find(truth_profile>nm+.5*contrast);
m.axial_span_um=0;
if ~isempty(above), m.axial_span_um=(above(end)-above(1)+1)*dx; end
m.axial_span_censored=~isempty(above) && (above(1)==1 || above(end)==sz(3));
m.true_axial_span_um=(true_above(end)-true_above(1)+1)*dx;
m.axial_span_error=(m.axial_span_um-m.true_axial_span_um)/m.true_axial_span_um;
segmented=n>nm+.5*contrast;
m.intersection_over_union=nnz(segmented & inside)/nnz(segmented | inside);
end
