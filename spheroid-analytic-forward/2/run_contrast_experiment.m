function run_contrast_experiment
%RUN_CONTRAST_EXPERIMENT Complete the three-index-contrast experiment here.
out=fileparts(mfilename('fullpath'));
addpath(fullfile(out,'..','1'));
parameter_sweep('contrast',out);
end
