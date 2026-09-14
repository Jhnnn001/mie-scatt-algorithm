function run_project(stage)
%RUN_PROJECT Reproduce the retained forward and inverse ODT experiments.
if nargin==0, stage='check'; end
stage=validatestring(stage,{'check','forward','fields','inverse'});
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'maxwell-solver'), ...
    fullfile(root,'born-rytov'),fullfile(root,'reconstruction'), ...
    fullfile(root,'experiments'), ...
    fullfile(root,'experiments','prolate-baseline'), ...
    fullfile(root,'experiments','size-sweep'), ...
    fullfile(root,'experiments','contrast-sweep'), ...
    fullfile(root,'experiments','field-analysis'));
switch stage
    case 'check'
        addpath(fullfile(root,'maxwell-solver','tests'), ...
            fullfile(root,'born-rytov','tests'),fullfile(root,'reconstruction','tests'));
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
        oblate_field_analysis('check');
    case 'forward'
        run_prolate_baseline;
        oblate_field_analysis('reference');
        run_size_experiment;
        run_contrast_experiment;
    case 'fields'
        oblate_field_analysis('xz');
        oblate_field_analysis('xz_validate');
    case 'inverse'
        base=struct('padding',8,'gp_iter',100,'tv_inner',200,'cg_max_iter',100);
        run_inverse_experiment('maxwell-direct-gp',base,{'direct','gp'});
        control=base; control.padding=4; control.source='matched';
        run_inverse_experiment('matched-control',control);
        tv=base; tv.methods={'tv'}; tv.tv_inner=1000;
        tv.tol=1e-5; tv.cg_tol=1e-5;
        % Archived penalty after matching the padding-4/8 physical normalization.
        tv.rho=0.003161909279393483;
        run_inverse_experiment('maxwell-tv',tv);
        check_inverse_sampling;
end
fprintf('PROJECT_STAGE_PASS %s\n',stage);
end
