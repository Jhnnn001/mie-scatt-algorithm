% Run the forward-model regressions from any working directory.
tests_dir = fileparts(mfilename('fullpath'));
forward_dir = fileparts(tests_dir);
addpath(forward_dir, tests_dir, fullfile(forward_dir,'..','spheroid-analytic-forward'));
runner_started = tic;
try
    test_born_rytov;
    fprintf('forward/tests/run_all: PASS (%.2f s)\n',toc(runner_started));
catch failure
    fprintf(2,'forward/tests/run_all: FAIL (%.2f s)\n',toc(runner_started));
    rethrow(failure);
end
