function result = demo_inverse(NA_det, model, opts)
%DEMO_INVERSE Water/cell spheroid -> existing forward -> three inverses.
% Example: result = demo_inverse(1,'rytov');
% Keep opts.NA_illum fixed when sweeping detector NA. No scalar NA is certified.
% Default length unit: um. No files are written by this function.
% Defaults use the selected low-contrast oblate. For exact-reference data,
% the full 512-square detector and a saved report, run_inverse_experiment.

if nargin < 1, NA_det = 1; end
if nargin < 2, model = 'rytov'; end
if nargin < 3, opts = struct; end
model = validatestring(model,{'born','rytov'});
defaults = struct('NA_illum',0.5,'n_angles',81,'grid_size',[128,128,80],'dx',0.1, ...
    'a',3,'b',5,'delta_scale',.01,'noise_std',0,'alpha',[], ...
    'max_iter',100,'outer_iter',5,'tol',1e-4,'plot',usejava('jvm'));
if ~isstruct(opts) || ~isscalar(opts) || ...
        any(~ismember(fieldnames(opts),fieldnames(defaults)))
    error('demo_inverse:InvalidOptions','Unknown or invalid demonstration options.');
end
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(opts,names{j}), opts.(names{j}) = defaults.(names{j}); end
end
lambda = 0.532; n_m = 1.335381534;
validateattributes(opts.delta_scale,{'numeric'},{'real','finite','scalar','positive'});
n_cell=n_m+opts.delta_scale*(1.365-n_m);
validateattributes(NA_det, {'numeric'}, {'real','finite','scalar','positive','<',n_m});
validateattributes(opts.NA_illum, {'numeric'}, ...
    {'real','finite','scalar','nonnegative','<=',NA_det});
validateattributes(opts.n_angles, {'numeric'}, {'finite','scalar','integer','positive'});
validateattributes(opts.dx, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(opts.a, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(opts.b, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(opts.noise_std, {'numeric'}, {'real','finite','scalar','nonnegative'});
validateattributes(opts.plot, {'logical'}, {'scalar'});
sz = opts.grid_size;
if isscalar(sz), sz = [sz,sz,sz]; end
validateattributes(sz, {'numeric'}, {'real','finite','integer','numel',3,'>=',4});
sz = double(sz(:).');
if any(mod(sz,2)) || opts.a>=(sz(3)/2-1)*opts.dx || ...
        opts.b>=(min(sz(1:2))/2-1)*opts.dx
    error('demo_inverse:InvalidGrid','Use an even grid with background surrounding the spheroid.');
end
inverse_dir = fileparts(mfilename('fullpath'));
addpath(inverse_dir,fullfile(inverse_dir,'..','forward'));
x = (-sz(1)/2:sz(1)/2-1)*opts.dx;
y = (-sz(2)/2:sz(2)/2-1)*opts.dx;
z = (-sz(3)/2:sz(3)/2-1)*opts.dx;
[X,Y,Z] = ndgrid(x,y,z);
inside = (X/opts.b).^2+(Y/opts.b).^2+(Z/opts.a).^2<1;
truth = n_m*ones(sz); truth(inside) = n_cell;
% Deterministic disk of illumination directions, including the axial wave.
radius = sqrt((0:opts.n_angles-1).'/max(opts.n_angles-1,1));
theta = asin(opts.NA_illum/n_m*radius);
phi = (0:opts.n_angles-1).'*(pi*(3-sqrt(5)));
z_det = (sz(3)/2+10)*opts.dx;
[born,rytov,u0,geom] = born_rytov_voxel(lambda,truth,n_m,opts.dx, ...
    theta,phi,NA_det,z_det);
if strcmp(model,'born'), U = u0+born; else, U = u0+rytov; end
if opts.noise_std>0
    U = U+opts.noise_std/sqrt(2)*(randn(size(U))+1i*randn(size(U)));
end
data = odt_prepare(U,u0,geom,model);
if isempty(opts.alpha)
    initial=data.operator.direct(data.y);
    opts.alpha=1e-3*max(abs(initial(:)));
end
solver_opts = struct('alpha',opts.alpha,'max_iter',opts.max_iter, ...
    'outer_iter',opts.outer_iter,'tol',opts.tol);
methods = {'direct','gp','tv'};
reconstructions = struct('method',{},'n',{},'chi',{},'info',{});
ri_error = zeros(3,1); mean_index = ri_error; axial_extent = ri_error;
sample_residual = ri_error; seconds = ri_error; converged = false(3,1);
for j = 1:3
    [n,chi,info] = odt_reconstruct(data,methods{j},solver_opts);
    reconstructions(j) = struct('method',methods{j},'n',n,'chi',chi,'info',info);
    ri_error(j) = norm(n(:)-truth(:))/norm(truth(:)-n_m);
    mean_index(j) = mean(n(inside)); % True support is used ONLY for evaluation.
    profile = squeeze(n(sz(1)/2+1,sz(2)/2+1,:));
    occupied = find(profile>n_m+0.5*(n_cell-n_m));
    if ~isempty(occupied), axial_extent(j) = (occupied(end)-occupied(1)+1)*opts.dx; end
    sample_residual(j) = info.relative_data_residual;
    seconds(j) = info.time; converged(j) = info.converged;
end
metrics = table(string(methods.'),ri_error,mean_index,axial_extent,sample_residual, ...
    seconds,converged,'VariableNames',{'method','relative_ri_error','mean_object_n', ...
    'axial_half_contrast_span_um','raw_sample_residual','seconds','converged'});
fprintf('%s: NA_det=%g, NA_illum=%g, %d angles, detector z=%g um\n', ...
    model,NA_det,opts.NA_illum,opts.n_angles,z_det);
fprintf('Physical post-pupil data; matched trilinear Ewald operator, %d samples.\n',numel(data.y));
disp(metrics);
result = struct('truth',truth,'reconstructions',reconstructions,'metrics',metrics, ...
    'data',data,'U',U,'u0',u0,'lambda',lambda,'n_m',n_m,'n_cell',n_cell, ...
    'NA_det',NA_det,'NA_illum',opts.NA_illum,'theta_inc',theta,'phi_inc',phi, ...
    'z_det',z_det,'x',x,'y',y,'z',z,'options',opts);
if opts.plot
    if ~usejava('jvm')
        error('demo_inverse:PlotUnavailable','Set opts.plot=false when MATLAB runs without a JVM.');
    end
    figure('Name',sprintf('%s ODT, detector NA %.3g',model,NA_det),'Color','w');
    volumes = [{truth},{reconstructions.n}];
    labels = [{'truth'},methods];
    for j = 1:4
        subplot(2,4,j);
        imagesc(x,y,volumes{j}(:,:,sz(3)/2+1).');
        axis image; axis xy; clim([n_m,n_cell]); colorbar;
        title(labels{j}); xlabel('x (um)'); ylabel('y (um)');
        subplot(2,4,j+4);
        imagesc(x,z,squeeze(volumes{j}(:,sz(2)/2+1,:)).');
        axis image; axis xy; clim([n_m,n_cell]); colorbar;
        xlabel('x (um)'); ylabel('z (um)');
    end
end
end
