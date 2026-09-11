% Run the inverse regressions without changing the caller's working directory.
test_dir = fileparts(mfilename('fullpath'));
addpath(test_dir, fileparts(test_dir), ...
    fullfile(test_dir,'..','..','forward'));
test_inverse;
test_matched_inverse;
test_exact_data;
