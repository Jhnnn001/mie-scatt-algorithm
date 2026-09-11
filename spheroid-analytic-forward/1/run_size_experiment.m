function run_size_experiment
%RUN_SIZE_EXPERIMENT Complete the three-size oblate experiment in this folder.
out=fileparts(mfilename('fullpath'));
addpath(out);
parameter_sweep('size',out);
end
