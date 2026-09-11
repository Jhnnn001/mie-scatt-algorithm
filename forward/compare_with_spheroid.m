%COMPARE_WITH_SPHEROID Scalar Born/Rytov versus co-polarized vector scattering.
% Run only after the complete spheroid-analytic-forward/tests/run_all.m is
% green; that external suite gate is separate from per-call info.validated.
% The vector first-order far field has co-polarized dipole factor
% 1-(rhat*E0').^2. It is one in the H plane and zero at E-plane 90 degrees;
% that zero need not persist in the exact finite-contrast vector solution.
% Near/interior scalar-vector differences need not vanish at weak contrast,
% especially for small x; those errors do not isolate Born/Rytov truncation.
% Detector-line log checks are practical branch/aliasing checks, not proof.
% Outputs remain in comparison_results, comparison_cases, comparison_check,
% comparison_figures, and comparison_elapsed; no binary files are written.

forward_dir = fileparts(mfilename('fullpath'));
addpath(forward_dir, fullfile(forward_dir,'..','spheroid-analytic-forward'));
comparison_started = tic;
n_m = 1;
lambda = 2*pi;
k = 2*pi*n_m/lambda;
shape_names = ["sphere", "prolate", "oblate"];
axis_factors = [1,1; 1,0.5; 0.5,1];
[contrast_grid, size_grid, shape_grid, incidence_grid] = ndgrid( ...
    [1e-3,1e-2,3e-2,1e-1], [1,5,20], 1:3, [0,pi/6]);
parameters = [contrast_grid(:),size_grid(:),shape_grid(:),incidence_grid(:)];
baseline_index = find(parameters(:,1)==1e-3 & parameters(:,2)==5 & ...
    parameters(:,3)==1 & parameters(:,4)==0);
parameters = parameters([baseline_index; ...
    setdiff((1:size(parameters,1)).',baseline_index,'stable')],:);

%% Evaluate the sole accuracy check first and retain its sweep case.
baseline = evaluate_case(lambda,n_m,1e-3,5,5,0,"sphere");
selection = baseline.indices.ring_H(baseline.theta<=pi/3);
low_errors = [relative_error(baseline.uB(selection),baseline.u_gt(selection)); ...
    relative_error(baseline.uR(selection),baseline.u_gt(selection))];
high = evaluate_fields(lambda,n_m,2e-3,5,5,0, ...
    baseline.points(baseline.indices.ring_H,:),baseline.E0);
selection = baseline.theta<=pi/3;
high_errors = [relative_error(high.uB(selection),high.u_gt(selection)); ...
    relative_error(high.uR(selection),high.u_gt(selection))];
error_ratios = high_errors./low_errors;
comparison_check = table(["Born";"Rytov"],low_errors,high_errors,error_ratios, ...
    'VariableNames',{'model','error_1e3','error_2e3','error_ratio'});
disp(comparison_check);
assert(all(low_errors<=0.02) && all(error_ratios>=1.5 & error_ratios<=2.5), ...
    'compare_with_spheroid:SphereAccuracy', ...
    'Sphere x=5 H-ring 0-60 deg must satisfy weak-contrast accuracy and scaling.');
fprintf('Sole sphere accuracy assertion: PASS (%.2f s)\n',toc(comparison_started));

%% Complete the 72-case sweep with one scalar and one vector call per case.
case_tables = cell(size(parameters,1),1);
comparison_cases = repmat(baseline,size(parameters,1),1);
for case_id = 1:size(parameters,1)
    p = parameters(case_id,:);
    shape_id = p(3);
    if case_id==1
        c = baseline;
    else
        axes_ab = (p(2)/k)*axis_factors(shape_id,:);
        c = evaluate_case(lambda,n_m,p(1),axes_ab(1),axes_ab(2),p(4), ...
            shape_names(shape_id));
    end
    [case_tables{case_id},c.phase_H,c.phase_E] = case_metrics(c,case_id);
    comparison_cases(case_id) = c;
    fprintf('%2d/72 %-7s x=%2g n_delta=%g zeta=%2g deg: %.2f s; H %s; E %s\n', ...
        case_id,c.shape,c.x,c.n_delta,c.zeta*180/pi,c.elapsed, ...
        c.phase_H.status,c.phase_E.status);
end
comparison_results = vertcat(case_tables{:});
disp(comparison_results);

%% Nine grouped figures: O1/O2/O3 for each representative shape.
comparison_figures = gobjects(3,3);
for shape_id = 1:3
    representative = find([comparison_cases.x]==20 & ...
        [comparison_cases.n_delta]==0.03 & [comparison_cases.zeta]==pi/6 & ...
        [comparison_cases.shape]==shape_names(shape_id),1);
    comparison_figures(shape_id,:) = plot_case(comparison_cases(representative));
end
comparison_elapsed = toc(comparison_started);
phase_rows = ismember(comparison_results.observable,["detector_H","detector_E"]);
fprintf('Comparison complete: %d cases, %d rows, %d figures, %.2f s.\n', ...
    numel(comparison_cases),height(comparison_results), ...
    numel(comparison_figures),comparison_elapsed);
fprintf('Detector-line phase metrics: %d accepted, %d skipped (two models per line).\n', ...
    nnz(phase_rows & comparison_results.phase_status=="accepted"), ...
    nnz(phase_rows & startsWith(comparison_results.phase_status,"skip:")));

function c = evaluate_case(lambda,n_m,n_delta,a,b,zeta,shape)
started = tic;
k = 2*pi*n_m/lambda;
R = max(a,b);
[~,~,khat,E0] = incident_plane_wave(k,n_m,zeta,0,[1,0],0,0,0);
e1 = E0;
e2 = cross(khat,e1);
s = linspace(-3*R,3*R,101).';
s_fine = linspace(-3*R,3*R,201).';
[s1,s2] = ndgrid(s,s);
detector = 2*R*khat+s1(:)*e1+s2(:)*e2;
line_H_fine = 2*R*khat+s_fine*e2;
line_E_fine = 2*R*khat+s_fine*e1;
theta = (0:180).'*pi/180;
rhat_H = cos(theta)*khat+sin(theta)*e2;
rhat_E = cos(theta)*khat+sin(theta)*e1;
axial_s = linspace(-3*R,3*R,601).';
axial_points = axial_s*khat;
rho_e = sqrt((axial_points(:,1).^2+axial_points(:,2).^2)/b^2 + ...
    axial_points(:,3).^2/a^2);
axial_keep = abs(rho_e-1)*min(a,b)>=0.02*R;
points = [detector;line_H_fine;line_E_fine;1e7*R*rhat_H; ...
    1e7*R*rhat_E;axial_points(axial_keep,:)];
ends = cumsum([0,size(detector,1),201,201,181,181,nnz(axial_keep)]);
indices = struct;
names = {'detector','fine_H','fine_E','ring_H','ring_E','axial'};
for j = 1:numel(names)
    indices.(names{j}) = (ends(j)+1:ends(j+1)).';
end
detector_indices = reshape(indices.detector,101,101);
indices.detector_H = detector_indices(51,:).';
indices.detector_E = detector_indices(:,51);
c = evaluate_fields(lambda,n_m,n_delta,a,b,zeta,points,E0);
c.shape = shape;
c.lambda = lambda;
c.n_m = n_m;
c.n_p = n_m*(1+n_delta);
c.n_delta = n_delta;
c.k = k;
c.x = k*R;
c.R = R;
c.a = a;
c.b = b;
c.zeta = zeta;
c.phi_inc = 0;
c.khat = khat;
c.E0 = E0;
c.e1 = e1;
c.e2 = e2;
c.Lk = 1/sqrt((khat(1)^2+khat(2)^2)/b^2+khat(3)^2/a^2);
c.DeltaPhi = 2*k*c.Lk*n_delta;
c.points = points;
c.indices = indices;
c.s = s;
c.s_fine = s_fine;
c.theta = theta;
c.rhat_H = rhat_H;
c.rhat_E = rhat_E;
c.dipole_E = 1-(rhat_E*E0').^2;
c.uB_E_dipole = c.uB(indices.ring_E).*c.dipole_E;
c.uR_E_dipole = c.uR(indices.ring_E).*c.dipole_E;
c.axial_s = axial_s;
c.axial_keep = axial_keep;
c.phase_H = struct;
c.phase_E = struct;
c.elapsed = toc(started);
end

function c = evaluate_fields(lambda,n_m,n_delta,a,b,zeta,points,E0)
n_p = n_m*(1+n_delta);
[uB,uR,u0,scalar_info] = born_rytov_spheroid(lambda,n_p,n_m,a,b,zeta,0, ...
    points(:,1),points(:,2),points(:,3));
opts = struct('field','scattered','tol',1e-6,'on_fail','error');
[Ex,Ey,Ez,vector_info] = spheroid_field(lambda,n_p,n_m,a,b,zeta,0,[1,0], ...
    points(:,1),points(:,2),points(:,3),opts);
if ~scalar_info.validated || ~vector_info.validated
    error('compare_with_spheroid:IneligibleCase', ...
        'Both solvers must validate every observation point before comparison.');
end
E = [Ex(:),Ey(:),Ez(:)];
if any(~isfinite(E),'all')
    error('compare_with_spheroid:NonfiniteGroundTruth', ...
        'The vector ground truth contains non-finite values.');
end
if isfield(vector_info,'sol')
    vector_info = rmfield(vector_info,'sol');
end
c = struct('uB',uB(:),'uR',uR(:),'u0',u0(:),'E',E,'u_gt',E*E0', ...
    'scalar_info',scalar_info,'vector_info',vector_info);
end

function [rows,phase_H,phase_E] = case_metrics(c,case_id)
names = ["detector";"detector_H";"detector_E";"ring_H";"ring_E";"axial";"ring_E_dipole"];
observable = repelem(names,2);
model = repmat(["Born";"Rytov"],numel(names),1);
relative_L2 = zeros(size(observable));
cross_fraction = zeros(size(observable));
phase_error = nan(size(observable));
phase_status = repmat("not_applicable",size(observable));
phase_H = line_phase(c,c.indices.detector_H,c.indices.fine_H);
phase_E = line_phase(c,c.indices.detector_E,c.indices.fine_E);
for j = 1:numel(names)
    pair = 2*j-1:2*j;
    if names(j)=="ring_E_dipole"
        index = c.indices.ring_E;
        models = [c.uB_E_dipole,c.uR_E_dipole];
    else
        index = c.indices.(names(j));
        models = [c.uB(index),c.uR(index)];
    end
    gt = c.u_gt(index);
    relative_L2(pair) = [relative_error(models(:,1),gt);relative_error(models(:,2),gt)];
    cross_fraction(pair) = norm_ratio(norm(c.E(index,:)-gt*c.E0,'fro'), ...
        norm(c.E(index,:),'fro'));
    if names(j)=="detector_H"
        phase_error(pair) = phase_H.errors;
        phase_status(pair) = phase_H.status;
    elseif names(j)=="detector_E"
        phase_error(pair) = phase_E.errors;
        phase_status(pair) = phase_E.status;
    end
end
n = numel(observable);
rows = table(repmat(case_id,n,1),repmat(c.shape,n,1),repmat(c.x,n,1), ...
    repmat(c.n_delta,n,1),repmat(c.a,n,1),repmat(c.b,n,1), ...
    repmat(c.zeta*180/pi,n,1),repmat(c.phi_inc,n,1), ...
    repmat(c.lambda,n,1),repmat(c.n_m,n,1),repmat(c.n_p,n,1), ...
    repmat(c.DeltaPhi,n,1),observable,model,relative_L2,cross_fraction, ...
    phase_error,phase_status,'VariableNames', ...
    {'case_id','shape','x','n_delta','a','b','zeta_deg','phi_inc', ...
    'lambda','n_m','n_p','DeltaPhi','observable','model','relative_L2', ...
    'cross_fraction','phase_error','phase_status'});
end

function phase = line_phase(c,coarse,fine)
endpoints = abs(c.u_gt(coarse([1,end]))./c.u0(coarse([1,end])));
reversed = endpoints(2)<endpoints(1);
if reversed
    coarse = flipud(coarse);
    fine = flipud(fine);
end
gt_ratio = c.u_gt(coarse)./c.u0(coarse);
gt_total = 1+gt_ratio;
fine_total = 1+c.u_gt(fine)./c.u0(fine);
coarse_phase = unwrap(angle(gt_total));
fine_phase = unwrap(angle(fine_total));
phase_difference = max(abs(coarse_phase-fine_phase(1:2:end)));
failed = strings(0,1);
if ~(abs(gt_ratio(1))<1e-2)
    failed(end+1) = "endpoint";
end
if ~all(abs(gt_total)>0.1) || ~all(abs(fine_total)>0.1)
    failed(end+1) = "total_amplitude";
end
if ~(phase_difference<=1e-3)
    failed(end+1) = "phase_refinement";
end
phase = struct('status',"accepted",'errors',[NaN;NaN], ...
    'reversed',reversed,'coarse_indices',coarse,'fine_indices',fine, ...
    'endpoint_ratio',abs(gt_ratio(1)), ...
    'minimum_total_amplitude',min(abs([gt_total;fine_total])), ...
    'refinement_difference',phase_difference,'phi_gt',[], ...
    'phi_Born',[],'phi_Rytov',[]);
if ~isempty(failed)
    phase.status = "skip: "+strjoin(failed,", ");
    return
end
phase.phi_gt = log(abs(gt_total))+1i*coarse_phase;
born_total = 1+c.uB(coarse)./c.u0(coarse);
phase.phi_Born = log(abs(born_total))+1i*unwrap(angle(born_total));
phase.phi_Rytov = c.uB(coarse)./c.u0(coarse);
phase.errors = [relative_error(phase.phi_Born,phase.phi_gt); ...
    relative_error(phase.phi_Rytov,phase.phi_gt)];
end

function value = relative_error(model,reference)
value = norm_ratio(norm(model-reference),norm(reference));
end

function value = norm_ratio(numerator,denominator)
if denominator==0
    if numerator==0
        value = 0;
    else
        value = Inf;
    end
else
    value = numerator/denominator;
end
end

function figures = plot_case(c)
label = sprintf('%s: x=%g, n_delta=%g, zeta=%g deg, DeltaPhi=%.3g rad', ...
    c.shape,c.x,c.n_delta,c.zeta*180/pi,c.DeltaPhi);
names = {'Vector co-pol.','Born','Rytov'};
fields = [c.u_gt,c.uB,c.uR];
figures = gobjects(1,3);
figures(1) = figure('Name',['O1 detector - ' label],'Color','w', ...
    'Position',[100,100,1200,720]);
layout = tiledlayout(2,6,'TileSpacing','compact');
title(layout,['O1 detector: ' label],'Interpreter','none');
color_max = max(abs(fields(c.indices.detector,:)),[],'all');
for j = 1:3
    nexttile([1,2]);
    imagesc(c.s/c.R,c.s/c.R,reshape(abs(fields(c.indices.detector,j)),101,101).');
    axis xy equal tight;
    if color_max>0
        clim([0,color_max]);
    end
    colorbar;
    title([names{j} ' |u_{sca}|']);
    xlabel('s_1/R (E)');
    ylabel('s_2/R (H)');
end
for plane = ["H","E"]
    nexttile([1,3]);
    plot(c.s/c.R,abs(fields(c.indices.("detector_"+plane),:)),'LineWidth',1);
    grid on;
    xlabel('s/R');
    ylabel('|u_{sca}|');
    title(plane+" central line");
    legend(names,'Location','best');
end

figures(2) = figure('Name',['O2 rings - ' label],'Color','w', ...
    'Position',[100,100,1200,460]);
layout = tiledlayout(1,2,'TileSpacing','compact');
title(layout,['O2 far rings: ' label],'Interpreter','none');
nexttile;
semilogy(c.theta*180/pi,abs(fields(c.indices.ring_H,:)),'LineWidth',1);
grid on;
xlabel('theta (deg)');
ylabel('|u_{sca}|');
title('H plane');
legend(names,'Location','best');
nexttile;
semilogy(c.theta*180/pi,abs([fields(c.indices.ring_E,:), ...
    c.uB_E_dipole,c.uR_E_dipole]),'LineWidth',1);
grid on;
xlabel('theta (deg)');
ylabel('|u_{sca}|');
title('E plane, including first-order dipole correction');
legend([names,{'Born x dipole','Rytov x dipole'}],'Location','best');

figures(3) = figure('Name',['O3 axis - ' label],'Color','w', ...
    'Position',[100,100,1000,750]);
layout = tiledlayout(3,1,'TileSpacing','compact');
title(layout,['O3 axial line: ' label],'Interpreter','none');
index = c.indices.axial;
total = c.u0(index)+fields(index,:);
kept = find(c.axial_keep);
breaks = [0;find(diff(kept)>1);numel(kept)];
accumulated_phase = zeros(size(total));
for segment = 1:numel(breaks)-1
    rows = breaks(segment)+1:breaks(segment+1);
    accumulated_phase(rows,:) = unwrap(angle(total(rows,:)./c.u0(index(rows))),[],1);
end
values = {real(total),imag(total),accumulated_phase};
ylabels = {'Re(u_{total})','Im(u_{total})','arg(u_{total}/u_0), unwrapped (rad)'};
for j = 1:3
    nexttile;
    plotted = nan(numel(c.axial_s),3);
    plotted(c.axial_keep,:) = values{j};
    plot(c.axial_s/c.R,plotted,'LineWidth',1);
    grid on;
    xlabel('axial distance/R');
    ylabel(ylabels{j});
    legend(names,'Location','best');
end
drawnow;
end
