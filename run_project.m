function run_project(stage)
%RUN_PROJECT Reproduce the retained forward and inverse ODT experiments.
if nargin==0, stage='check'; end
stage=validatestring(stage,{'check','forward','fields','inverse'});
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'spheroid-analytic-forward'), ...
    fullfile(root,'spheroid-analytic-forward','1'), ...
    fullfile(root,'spheroid-analytic-forward','2'), ...
    fullfile(root,'forward'),fullfile(root,'forward','exp'), ...
    fullfile(root,'inverse'));
switch stage
    case 'check'
        addpath(fullfile(root,'spheroid-analytic-forward','tests'), ...
            fullfile(root,'forward','tests'),fullfile(root,'inverse','tests'));
        test_gauss_legendre;
        test_incident_plane_wave;
        test_sph_coords;
        test_special_functions;
        test_coords_vecwave;
        test_born_rytov;
        test_born_rytov_voxel;
        test_inverse;
        test_matched_inverse;
        test_exact_data;
        oblate_na_sweep('check');
    case 'forward'
        addpath(fullfile(root,'forward','tests'));
        analyze_rytov_spheroid;
        oblate_na_sweep('reference');
        run_size_experiment;
        run_contrast_experiment;
    case 'fields'
        oblate_na_sweep('xz');
        oblate_na_sweep('xz_validate');
    case 'inverse'
        base=struct('padding',8,'gp_iter',100,'tv_inner',200,'cg_max_iter',100);
        run_inverse_experiment('padding8',base,{'direct','gp'});
        control=base; control.padding=4; control.source='matched';
        run_inverse_experiment('linear_control',control);
        tv=base; tv.methods={'tv'}; tv.tv_inner=1000;
        tv.tol=1e-5; tv.cg_tol=1e-5;
        % Archived penalty after matching the padding-4/8 physical normalization.
        tv.rho=0.003161909279393483;
        run_inverse_experiment('tv_strict',tv);
        check_inverse_sampling;
end
fprintf('PROJECT_STAGE_PASS %s\n',stage);
end
