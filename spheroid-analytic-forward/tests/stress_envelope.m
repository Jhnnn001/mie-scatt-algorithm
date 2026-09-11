%STRESS_ENVELOPE Deterministic coverage report for the tested envelope.
% Set stress_output_dir before running to override the default report path.

test_directory = fileparts(mfilename('fullpath'));
solver_directory = fileparts(test_directory);
repository_directory = fileparts(solver_directory);
addpath(solver_directory);

if ~exist('stress_output_dir', 'var') || isempty(stress_output_dir)
    stress_output_dir = fullfile(repository_directory, '..', 'tmp', ...
        'spheroid-analytic-forward', 'stress');
end
if isstring(stress_output_dir) && isscalar(stress_output_dir)
    stress_output_dir = char(stress_output_dir);
end
if ~(ischar(stress_output_dir) && isrow(stress_output_dir) && ...
        ~isempty(stress_output_dir))
    error('stress_envelope:InvalidOutputDirectory', ...
        'stress_output_dir must be a nonempty character row or string scalar.');
end
if exist(stress_output_dir, 'dir') ~= 7
    [created, message] = mkdir(stress_output_dir);
    if ~created
        error('stress_envelope:OutputDirectory', '%s', message);
    end
end

stress_cases = build_cases();
report_path = fullfile(stress_output_dir, 'stress_report.mat');
csv_path = fullfile(stress_output_dir, 'stress_report.csv');

if exist(report_path, 'file') == 2
    loaded_report = load(report_path, 'stress_report');
    if ~isfield(loaded_report, 'stress_report') || ...
            ~compatible_report(loaded_report.stress_report, stress_cases)
        error('stress_envelope:IncompatibleReport', ...
            'The existing MAT report describes a different campaign.');
    end
    stress_report = loaded_report.stress_report;
else
    stress_report = struct('schema_version', 1, ...
        'description', ['Section 6 corners, 20 deterministic ', ...
        'Latin-hypercube samples, and the x=30:0.01:35 resonance sweep'], ...
        'started_at', timestamp(), 'updated_at', '', ...
        'total_case_count', numel(stress_cases), 'completed_count', 0, ...
        'minimum_corner_x_ext', 0.02, 'cases', stress_cases, ...
        'rows', repmat(empty_row(), 0, 1));
end
completed_ids = {stress_report.rows.case_id};
if numel(unique(completed_ids)) ~= numel(completed_ids)
    error('stress_envelope:DuplicateCase', ...
        'The existing report contains duplicate case identifiers.');
end
if ~all(ismember(completed_ids, {stress_cases.case_id}))
    error('stress_envelope:UnknownCase', ...
        'The existing report contains an unknown case identifier.');
end
% The MAT file is authoritative; regenerate a stale CSV before resuming.
save_report(stress_report, report_path, csv_path);
fprintf('stress_envelope: %d/%d cases already recorded in %s\n', ...
    numel(completed_ids), numel(stress_cases), stress_output_dir);
completed_cases = ismember({stress_cases.case_id}, completed_ids);

solver_options = struct('tol', 1e-6, 'Lmax', 400, 'Mmax', 300, ...
    'on_fail', 'warn', 'verbose', false);
for case_index = 1:numel(stress_cases)
    case_spec = stress_cases(case_index);
    if completed_cases(case_index)
        continue
    end

    row = empty_row();
    row = copy_parameters(row, case_spec, case_index, solver_options.tol);
    started = tic;
    lastwarn('');
    try
        sol = spheroid_solve(case_spec.lambda, case_spec.n_p, ...
            case_spec.n_m, case_spec.a, case_spec.b, ...
            case_spec.theta_inc, case_spec.phi_inc, case_spec.pol, ...
            solver_options);
        row = copy_result(row, sol.info);
        row.status = 'returned';
    catch exception
        if strcmp(exception.identifier, 'MATLAB:OperationTerminated')
            rethrow(exception)
        end
        row.status = 'error';
        row.error_identifier = exception.identifier;
        row.error_message = exception.message;
        row.error_detail = struct('identifier', exception.identifier, ...
            'message', exception.message, 'stack', exception.stack);
    end
    row.elapsed_seconds = toc(started);
    [row.warning_message, row.warning_identifier] = lastwarn;
    row.completed_at = timestamp();
    stress_report.rows(end + 1, 1) = row;
    stress_report.completed_count = numel(stress_report.rows);
    stress_report.updated_at = row.completed_at;
    save_report(stress_report, report_path, csv_path);
    completed_cases(case_index) = true;

    fprintf('[%d/%d] %-13s %-16s %8.2f s  validated=%d\n', ...
        case_index, numel(stress_cases), case_spec.family, ...
        row.status, row.elapsed_seconds, row.validated);
end

fprintf('stress_envelope: report saved to %s and %s\n', ...
    report_path, csv_path);

function cases = build_cases()
template = case_template();
sample_count = 20;
resonance_x = (30:0.01:35).';
cases = repmat(template, 16 + sample_count + numel(resonance_x), 1);
next_case = 0;
corner_index = 0;
for size_level = 1:2
    % The plan's oblate xi0 >= 0.26 note rounds the exact a/b = 1/4
    % endpoint, for which xi0 = 1/sqrt(15).
    for aspect = [1/4, 4]
        for n_real = [0.5, 2]
            for loss_ratio = [0, 0.1]
                corner_index = corner_index + 1;
                n_relative = n_real*(1 + 1i*loss_ratio);
                if size_level == 1
                    x_ext = 0.02;
                else
                    x_ext = 250/max(1, abs(n_relative));
                end
                next_case = next_case + 1;
                cases(next_case) = make_case( ...
                    sprintf('corner_%03d', corner_index), 'corner', ...
                    corner_index, x_ext, aspect, n_relative);
            end
        end
    end
end

stream = RandStream('mt19937ar', 'Seed', 20260912);
unit_samples = zeros(sample_count, 4);
for dimension = 1:size(unit_samples, 2)
    permutation = randperm(stream, sample_count).';
    unit_samples(:, dimension) = ...
        (permutation - rand(stream, sample_count, 1))/sample_count;
end
for sample_index = 1:sample_count
    aspect = exp(log(1/4) + unit_samples(sample_index, 2)*log(16));
    n_real = 0.5 + 1.5*unit_samples(sample_index, 3);
    loss_ratio = 0.1*unit_samples(sample_index, 4);
    n_relative = n_real*(1 + 1i*loss_ratio);
    x_limit = 250/max(1, abs(n_relative));
    x_ext = exp(log(0.02) + unit_samples(sample_index, 1)* ...
        (log(x_limit) - log(0.02)));
    next_case = next_case + 1;
    cases(next_case) = make_case( ...
        sprintf('lhs_%03d', sample_index), 'latin_hypercube', ...
        sample_index, x_ext, aspect, n_relative);
end

for sweep_index = 1:numel(resonance_x)
    next_case = next_case + 1;
    cases(next_case) = make_case( ...
        sprintf('resonance_%03d', sweep_index), 'resonance', ...
        sweep_index, resonance_x(sweep_index), 1.2, 1.5);
end
assert(next_case == numel(cases));
end

function value = case_template()
value = struct('case_id', '', 'family', '', 'family_index', NaN, ...
    'lambda', NaN, 'n_p', complex(NaN), 'n_m', NaN, 'a', NaN, ...
    'b', NaN, 'theta_inc', NaN, 'phi_inc', NaN, ...
    'pol', complex([NaN, NaN]), 'x_ext', NaN, 'x_int', NaN, ...
    'aspect', NaN, 'loss_ratio', NaN);
end

function value = make_case(case_id, family, family_index, x_ext, ...
    aspect, n_relative)
value = case_template();
value.case_id = case_id;
value.family = family;
value.family_index = family_index;
value.lambda = 2*pi;
value.n_p = n_relative;
value.n_m = 1;
if aspect > 1
    value.a = x_ext;
    value.b = x_ext/aspect;
else
    value.a = aspect*x_ext;
    value.b = x_ext;
end
value.theta_inc = pi/6;
value.phi_inc = 0;
value.pol = [0, 1];
value.x_ext = x_ext;
value.x_int = abs(n_relative)*x_ext;
value.aspect = aspect;
value.loss_ratio = imag(n_relative)/real(n_relative);
end

function row = empty_row()
row = struct('case_index', NaN, 'case_id', '', 'family', '', ...
    'family_index', NaN, 'lambda', NaN, 'n_m', NaN, ...
    'n_p_real', NaN, 'n_p_imag', NaN, 'n_relative_real', NaN, ...
    'n_relative_imag', NaN, 'loss_ratio', NaN, 'aspect', NaN, ...
    'a', NaN, 'b', NaN, 'x_ext', NaN, 'x_int', NaN, ...
    'theta_inc', NaN, 'phi_inc', NaN, 'pol_1_real', NaN, ...
    'pol_1_imag', NaN, 'pol_2_real', NaN, 'pol_2_imag', NaN, ...
    'tolerance', NaN, 'status', '', 'elapsed_seconds', NaN, ...
    'type', '', 'L', NaN, 'M', NaN, 'Qtheta', NaN, 'Qphi', NaN, ...
    'rmax', NaN, 'xi_switch', NaN, 'xi_start', NaN, ...
    'residual', NaN, 'rcond_min', NaN, 'rcond_min_active', NaN, ...
    'rounds', NaN, 'inside_envelope', false, ...
    'all_checks_passed', false, 'validated', false, ...
    'failure_reason', '', 'warning_identifier', '', ...
    'warning_message', '', 'error_identifier', '', 'error_message', '', ...
    'completed_at', '', 'checks', struct(), ...
    'error_detail', struct('identifier', '', 'message', '', ...
    'stack', []));
names = check_names();
for index = 1:numel(names)
    row.([names{index}, '_computed']) = false;
    row.([names{index}, '_value']) = NaN;
    row.([names{index}, '_passed']) = false;
end
end

function row = copy_parameters(row, case_spec, case_index, tolerance)
row.case_index = case_index;
row.case_id = case_spec.case_id;
row.family = case_spec.family;
row.family_index = case_spec.family_index;
row.lambda = case_spec.lambda;
row.n_m = case_spec.n_m;
row.n_p_real = real(case_spec.n_p);
row.n_p_imag = imag(case_spec.n_p);
row.n_relative_real = real(case_spec.n_p/case_spec.n_m);
row.n_relative_imag = imag(case_spec.n_p/case_spec.n_m);
row.loss_ratio = case_spec.loss_ratio;
row.aspect = case_spec.aspect;
row.a = case_spec.a;
row.b = case_spec.b;
row.x_ext = case_spec.x_ext;
row.x_int = case_spec.x_int;
row.theta_inc = case_spec.theta_inc;
row.phi_inc = case_spec.phi_inc;
row.pol_1_real = real(case_spec.pol(1));
row.pol_1_imag = imag(case_spec.pol(1));
row.pol_2_real = real(case_spec.pol(2));
row.pol_2_imag = imag(case_spec.pol(2));
row.tolerance = tolerance;
end

function row = copy_result(row, info)
row.type = info.type;
row.L = info.L;
row.M = info.M;
row.Qtheta = info.Qtheta;
row.Qphi = info.Qphi;
row.rmax = info.rmax;
row.xi_switch = info.xi_switch;
row.xi_start = info.xi_start;
row.residual = info.residual;
row.rcond_min = info.rcond_min;
row.rcond_min_active = info.rcond_min_active;
row.rounds = info.rounds;
row.inside_envelope = info.inside_envelope;
row.validated = info.validated;
row.failure_reason = info.failure_reason;
row.checks = info.checks;

names = check_names();
row.all_checks_passed = true;
for index = 1:numel(names)
    name = names{index};
    check = info.checks.(name);
    row.([name, '_computed']) = check.computed;
    row.([name, '_value']) = check.value;
    row.([name, '_passed']) = check.passed;
    row.all_checks_passed = row.all_checks_passed && ...
        check.computed && check.passed;
end
end

function names = check_names()
names = {'incident_tail', 'qphi_independence', 'modal_tail', ...
    'surface_residual', 'qtheta_independence', 'cutoff_stability', ...
    'special_functions'};
end

function value = timestamp()
value = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmss'));
end

function valid = compatible_report(report, stress_cases)
required = {'schema_version', 'cases', 'rows', 'total_case_count', ...
    'completed_count'};
valid = isstruct(report) && isscalar(report) && all(isfield(report, required));
if ~valid
    return
end
valid = isequal(report.schema_version, 1) && ...
    isequaln(report.cases, stress_cases) && isstruct(report.rows) && ...
    isfield(report.rows, 'case_id') && ...
    all(arrayfun(@(row) ischar(row.case_id) && isrow(row.case_id), ...
    report.rows)) && ...
    isequal(report.total_case_count, numel(stress_cases)) && ...
    isequal(report.completed_count, numel(report.rows));
end

function save_report(stress_report, report_path, csv_path)
temporary_mat = fullfile(fileparts(report_path), ...
    'stress_report.partial.mat');
save(temporary_mat, 'stress_report', '-v7');
if ~isempty(stress_report.rows)
    flat_rows = rmfield(stress_report.rows, {'checks', 'error_detail'});
    report_table = struct2table(flat_rows(:), 'AsArray', true);
    temporary_csv = fullfile(fileparts(csv_path), ...
        'stress_report.partial.csv');
    writetable(report_table, temporary_csv);
end
[moved, message] = movefile(temporary_mat, report_path, 'f');
if ~moved
    error('stress_envelope:SaveReport', '%s', message);
end
if ~isempty(stress_report.rows)
    [moved, message] = movefile(temporary_csv, csv_path, 'f');
    if ~moved
        error('stress_envelope:SaveReport', '%s', message);
    end
elseif exist(csv_path, 'file') == 2
    delete(csv_path);
end
end
