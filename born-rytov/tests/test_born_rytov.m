% Independent scalar Born references and Rytov regression checks.
test_started = tic;

% T4 catches missing Rytov mapping, cancellation, and non-finite acceptance.
points = [0,0,0; 0.2,0.3,0.1; 1.2,0.4,-0.2; -2,1,3];
[uB, uR, ~, info] = born_rytov_spheroid(2*pi, 1, 1, 1, 1, ...
    pi/6, pi/4, points(:,1), points(:,2), points(:,3));
assert(all(uB(:)==0) && all(uR(:)==0) && info.validated, ...
    'T4 zero contrast must return exact zero Born and Rytov fields.');
[uB, uR, u0, info] = born_rytov_spheroid(2*pi/5, 1.1, 1, 1, 1, ...
    pi/6, pi/4, points(:,1), points(:,2), points(:,3));
phi_R = uB./u0;
selected = abs(phi_R)>1e-3;
assert(any(selected), 'T4 finite-phase identity path not exercised.');
total = u0.*exp(phi_R);
assert(info.validated && all(abs(u0(selected)+uR(selected)-total(selected)) ...
    <= 1e-14*abs(total(selected))), 'T4 Rytov identity failed.');
far_angle = 5e-3;
far_point = 1e9*[sin(far_angle),0,cos(far_angle)];
[far_uB, far_uR, far_u0, far_info] = born_rytov_spheroid(2*pi/100,200,1,1,1, ...
    0,0,far_point(1),far_point(2),far_point(3));
far_expected = far_u0.*expm1(far_uB./far_u0);
far_identity_error = abs(far_uR-far_expected)/max(abs(far_uR),realmin);
assert(abs(far_uB/far_u0)>1e-2 && abs(far_uB/far_u0)<1 && ...
    far_info.validated && far_identity_error<=1e-12, ...
    'T4 far Rytov output is inconsistent with the returned Born and incident fields.');
rng(20260917, 'twister');
directions = randn(20,3);
directions = directions./vecnorm(directions,2,2);
points = directions.*(1.2+1.8*rand(20,1));
[uB, uR, u0, info] = born_rytov_spheroid(2*pi, 1+1e-9, 1, 1, 1, ...
    0, 0, points(:,1), points(:,2), points(:,3));
phi_R = uB./u0;
small_reference = phi_R+phi_R.^2/2;
assert(any(abs(exp(phi_R)-1-small_reference)>1e-8*abs(phi_R)), ...
    'T4 exp(phi)-1 cancellation path not exercised.');
assert(info.validated && all(abs(uR./u0-small_reference)<=1e-8*abs(phi_R)), ...
    'T4 small-phase Rytov lost accuracy.');
directions = randn(100,3);
directions = directions./vecnorm(directions,2,2);
points = directions.*(1.2+1.8*rand(100,1));
[uB, uR, ~, info] = born_rytov_spheroid(2*pi/5, 1+1e-4, 1, 1, 1, ...
    0, 0, points(:,1), points(:,2), points(:,3));
assert(info.validated && all(abs(uR-uB)<=1e-3*abs(uB)), ...
    'T4 weak-contrast Rytov disagrees with Born.');

overflow_point = [0,0,1.05];
overflow_args = {2*pi/30, 1+1e6, 1, 1, 1, 0, 0, ...
    overflow_point(1), overflow_point(2), overflow_point(3)};
assert_error(@() born_rytov_spheroid(overflow_args{:}), ...
    'born_rytov_spheroid:ConvergenceFailure');
warning_identifier = 'born_rytov_spheroid:ConvergenceFailure';
warning_state = warning('query',warning_identifier);
warning('on',warning_identifier);
warning_cleanup = onCleanup(@() warning(warning_state.state,warning_identifier));
lastwarn('');
overflow_warning_text = evalc('[overflow_uB, overflow_uR, overflow_u0, overflow_info] = born_rytov_spheroid(overflow_args{:}, struct(''on_fail'',''warn''));');
[~, warning_id] = lastwarn;
clear warning_cleanup
assert(isscalar(strfind(overflow_warning_text,'Warning:')), ...
    'T4 overflow warn policy must emit exactly one warning.');
assert(real(overflow_uB/overflow_u0)>log(realmax), ...
    'T4 exponential overflow path not exercised.');
assert(strcmp(warning_id,'born_rytov_spheroid:ConvergenceFailure') && ...
    isfinite(overflow_uB) && ~isfinite(overflow_uR) && ...
    overflow_info.unconverged && ~overflow_info.validated, ...
    'T4 non-finite Rytov did not follow the warn policy.');
fprintf('T4 zero contrast, identity, small phase, weak limit, and overflow PASS\n');

% This pair passes Born alone, but exponential amplification needs refinement.
amp_k = 30;
amp_ratio = sqrt(101);
amp_point = [0,0,1.05];
amp_order = 48;
amp_check = ceil(1.5*amp_order);
amp_tol = 1e-8;
amp_f = amp_k^2*(amp_ratio^2-1);
amp_u0 = exp(1i*amp_k*amp_point(3));
amp_coarse = naive_born_surface(1,amp_k,amp_f,amp_point,amp_order);
amp_fine = naive_born_surface(1,amp_k,amp_f,amp_point,amp_check);
amp_natural_scale = abs(amp_f)/(3*norm(amp_point));
amp_scale = max(abs(amp_fine),amp_natural_scale);
amp_r_coarse = amp_u0*expm1(amp_coarse/amp_u0);
amp_r_fine = amp_u0*expm1(amp_fine/amp_u0);
amp_r_scale = max(abs(amp_r_fine),amp_natural_scale);
assert(abs(amp_fine-amp_coarse)<=amp_tol*amp_scale && ...
    abs(amp_r_fine-amp_r_coarse)>amp_tol*amp_r_scale && ...
    isfinite(amp_r_coarse) && isfinite(amp_r_fine), ...
    'T4 finite Rytov amplification path not exercised.');
[amp_reference_short, amp_reference] = born_sphere_series(1,amp_k,amp_f,[0,0,1], ...
    amp_point,ceil(amp_k+4.05*amp_k^(1/3)+20),false);
amp_r_reference = amp_u0*expm1(amp_reference/amp_u0);
assert(all(isfinite(amp_r_reference(:))) && real(amp_reference/amp_u0)<log(realmax), ...
    'T4 amplification reference must remain finite.');
assert(abs(amp_r_reference-amp_u0*expm1(amp_reference_short/amp_u0)) ...
    <= 1e-10*abs(amp_r_reference), 'T4 amplified sphere reference tail did not converge.');
[~, amp_uR, ~, amp_info] = born_rytov_spheroid(2*pi/amp_k,amp_ratio,1,1,1, ...
    0,0,amp_point(1),amp_point(2),amp_point(3), ...
    struct('tol',amp_tol,'N_theta',amp_order));
amp_error = abs(amp_uR-amp_r_reference)/max(abs(amp_r_reference),amp_natural_scale);
assert(amp_info.validated && amp_info.N_theta>amp_check, ...
    'T4 Born-only acceptance missed Rytov refinement.');
assert(amp_error<=amp_tol, 'T4 amplified Rytov disagrees with the sphere reference.');
fprintf('T4 amplification: pair %d/%d, accepted %d, relative error %.3e PASS\n', ...
    amp_order,amp_check,amp_info.N_theta,amp_error);

% T5 catches a missing adaptive loop and the final cap-pair failure rule.
failure_opts = struct('N_theta', 76, 'N_max', 200);
failure_call = @() born_rytov_spheroid(2*pi/30, 1.1, 1, 1, 1, 0, 0, ...
    1+1e-4, 0, 0, failure_opts);
failure = assert_error(failure_call, 'born_rytov_spheroid:ConvergenceFailure');
assert(contains(failure.message, '1 point') && contains(failure.message, '0.0001'), ...
    'T5 convergence error must report the failed count and surface distance.');
failure_opts.on_fail = 'warn';
warning_state = warning('query',warning_identifier);
warning('on',warning_identifier);
warning_cleanup = onCleanup(@() warning(warning_state.state,warning_identifier));
lastwarn('');
warning_text = evalc('[failed_uB, ~, ~, failed_info] = born_rytov_spheroid(2*pi/30, 1.1, 1, 1, 1, 0, 0, 1+1e-4, 0, 0, failure_opts);');
[~, warning_id] = lastwarn;
lastwarn('');
mixed_warning_text = evalc('[mixed_uB, mixed_uR, ~, mixed_info] = born_rytov_spheroid(2*pi/30, 1.1, 1, 1, 1, 0, 0, [1+1e-4,3], [0,0], [0,0], failure_opts);');
[~, mixed_warning_id] = lastwarn;
clear warning_cleanup
assert(strcmp(warning_id, 'born_rytov_spheroid:ConvergenceFailure'), ...
    'T5 warn policy did not emit the convergence warning.');
assert(isscalar(strfind(warning_text, 'Warning:')), ...
    'T5 warn policy must emit exactly one warning.');
assert(strcmp(mixed_warning_id,warning_identifier) && ...
    isscalar(strfind(mixed_warning_text,'Warning:')), ...
    'T5 mixed warn policy must emit exactly one convergence warning.');
assert(failed_info.unconverged && ~failed_info.validated, 'T5 failure path not exercised.');
assert(failed_info.N_theta == 200, 'T5 did not evaluate the final cap pair.');
coarse = naive_born_surface(1, 30, 30^2*(1.1^2-1), [1+1e-4, 0, 0], 133);
fine = naive_born_surface(1, 30, 30^2*(1.1^2-1), [1+1e-4, 0, 0], 200);
natural_scale = abs(failed_info.f)*failed_info.V/(4*pi*(1+1e-4));
scale = max(abs(fine), natural_scale);
assert(abs(failed_uB-fine) <= 1e-12*scale, 'T5 failed value is not the last finer value.');
rytov_coarse = expm1(coarse);
rytov_fine = expm1(fine);
expected_error = max(abs(fine-coarse)/scale, ...
    abs(rytov_fine-rytov_coarse)/max(abs(rytov_fine),natural_scale));
assert(abs(failed_info.err-expected_error) <= 1e-12, ...
    'T5 info.err is not the normalized final pair difference.');
assert(failed_info.err > 1e-8, 'T5 failure path not exercised.');
adjacent_opts = struct('N_theta',76,'N_max',115,'on_fail','warn');
evalc('[adjacent_uB, ~, ~, adjacent_info] = born_rytov_spheroid(2*pi/30, 1.1, 1, 1, 1, 0, 0, 1+1e-4, 0, 0, adjacent_opts);');
adjacent_coarse = naive_born_surface(1,30,30^2*(1.1^2-1),[1+1e-4,0,0],76);
adjacent_fine = naive_born_surface(1,30,30^2*(1.1^2-1),[1+1e-4,0,0],115);
adjacent_scale = max(abs(adjacent_fine),natural_scale);
adjacent_rytov_coarse = expm1(adjacent_coarse);
adjacent_rytov_fine = expm1(adjacent_fine);
adjacent_expected_error = max(abs(adjacent_fine-adjacent_coarse)/adjacent_scale, ...
    abs(adjacent_rytov_fine-adjacent_rytov_coarse)/ ...
    max(abs(adjacent_rytov_fine),natural_scale));
assert(adjacent_info.unconverged, 'T5 adjacent-cap failure path not exercised.');
assert(adjacent_info.N_theta==115, 'T5 adjacent cap did not evaluate N_max.');
assert(abs(adjacent_uB-adjacent_fine)<=1e-12*adjacent_scale, ...
    'T5 adjacent cap did not return the N_max value.');
assert(abs(adjacent_info.err-adjacent_expected_error)<=1e-12, ...
    'T5 adjacent-cap info.err does not describe the N_theta/N_max pair.');
[far_only_uB, far_only_uR] = born_rytov_spheroid(2*pi/30,1.1,1,1,1,0,0,3,0,0);
assert(isequal(mixed_info.unconverged,[true,false]) && ~mixed_info.validated && ...
    all(isfinite(mixed_uB)) && all(isfinite(mixed_uR)) && ...
    abs(mixed_uB(2)-far_only_uB)<=1e-12*max(abs(far_only_uB),realmin) && ...
    abs(mixed_uR(2)-far_only_uR)<=1e-12*max(abs(far_only_uR),realmin), ...
    'T5 mixed warn result did not preserve the independently converged point.');
for policy = {'error', 'warn'}
    assert_error(@() born_rytov_spheroid(2*pi, 1.01, 1, 1, 1, 0, 0, ...
        1+1e-10, 0, 0, struct('on_fail', policy{1})), ...
        'born_rytov_spheroid:SurfacePoint');
end
assert_error(@() born_rytov_spheroid(2*pi, 1.01, 1, 1, 1, 0, 0, ...
    2, 0, 0, struct('N_theta', 76, 'N_max', 113)), ...
    'born_rytov_spheroid:InvalidQuadratureLimit');
assert_error(@() born_rytov_spheroid(2*pi, realmax, 1, 1, 1, 0, 0, ...
    2, 0, 0), 'born_rytov_spheroid:InvalidDerivedParameters');
assert_error(@() born_rytov_spheroid(2*pi,1.01,1,1,1,0,0,2,0,0, ...
    struct('tol',99*eps)), 'born_rytov_spheroid:InvalidTolerance');
[~, ~, ~, tolerance_floor_info] = born_rytov_spheroid(2*pi,1,1,1,1,0,0, ...
    2,0,0,struct('tol',100*eps));
assert(tolerance_floor_info.validated, ...
    'T5 must accept the documented tolerance floor.');
fprintf('T5 failure policy, final-pair metadata, and input rejection PASS\n');

% T1 catches the interior jump and adaptive Born accuracy independently.
rng(20260914, 'twister');
theta_inc = 40*pi/180;
phi_inc = 20*pi/180;
khat = [sin(theta_inc)*cos(phi_inc), sin(theta_inc)*sin(phi_inc), cos(theta_inc)];
sphere_cases = [0.5, 1.05; 5, 1.01; 30, 1.02+0.01i];
for case_index = 1:size(sphere_cases,1)
    k = real(sphere_cases(case_index,1));
    n_p = sphere_cases(case_index,2);
    f = k^2*(n_p^2-1);
    directions = randn(100, 3);
    directions = directions ./ vecnorm(directions, 2, 2);
    radii = [0.9*rand(50,1).^(1/3); 1.1+8.9*rand(50,1)];
    points = [directions.*radii; 0,0,0; zeros(5,2), [-3;-0.8;0.3;1.2;5]];
    lmax = ceil(k+4.05*k^(1/3)+20);
    [reference, reference_extended] = born_sphere_series(1, k, f, khat, points, lmax, false);
    reference_scale = max(abs(reference_extended));
    tail = max(abs(reference_extended-reference))/reference_scale;
    assert(tail < 1e-12, 'T1 sphere reference tail did not converge.');
    [uB, ~, ~, info] = born_rytov_spheroid(2*pi/k, n_p, 1, 1, 1, ...
        theta_inc, phi_inc, points(:,1), points(:,2), points(:,3));
    sphere_error = max(abs(uB-reference_extended))/reference_scale;
    assert(sphere_error <= 1e-8, 'T1 Born disagrees with the independent sphere series.');
    assert(info.validated, 'T1 Born quadrature was not validated.');
    fprintf('T1 x=%g: series tail %.3e, Born relative error %.3e PASS\n', k, tail, sphere_error);
end

% Exact finite Hankel sums retain all radial corrections at distant points.
rng(20260915, 'twister');
far_radii = [1e4,1e7,1e9];
for x = [5, 100]
    lambda = 2*pi/x;
    k = 2*pi/lambda;
    directions = [0,0,1; randn(19,3)];
    directions = directions ./ vecnorm(directions,2,2);
    points = [1e4*directions; 1e7*directions; 1e9*directions];
    rho = sqrt(sum(points.^2,2));
    f = k^2*(1.01^2-1);
    lmax = ceil(k+4.05*k^(1/3)+20);
    [reference, reference_extended] = born_sphere_series(1, k, f, [0,0,1], points, lmax, true);
    [uB, ~, ~, info] = born_rytov_spheroid(lambda, 1.01, 1, 1, 1, ...
        0, 0, points(:,1), points(:,2), points(:,3));
    dephased = uB.*exp(-1i*k*rho);
    for radius_index = 1:3
        index = (1:20)+20*(radius_index-1);
        reference_scale = max(abs(reference_extended(index)));
        tail = max(abs(reference(index)-reference_extended(index)))/reference_scale;
        far_error = max(abs(dephased(index)-reference_extended(index)))/reference_scale;
        assert(tail < 1e-12, 'T1 far sphere reference tail did not converge.');
        assert(far_error <= 1e-8, 'T1 dephased Born disagrees with the exact Hankel sum.');
        fprintf('T1 x=%g, r/a=%g: dephased relative error %.3e PASS\n', ...
            x, far_radii(radius_index), far_error);
    end
    assert(info.validated, 'T1 distant Born quadrature was not validated.');
    if x == 100
        naive = naive_born_surface(1, k, f, [0,0,1e9], 324)*exp(-1i*k*1e9);
        naive_error = abs(naive-reference_extended(41))/abs(reference_extended(41));
        assert(naive_error > 1e-8, 'T1 naive kernel failure path not exercised.');
        fprintf('T1 naive phase kernel relative error %.3e exceeds 1e-8 as expected\n', naive_error);
    end
end

% T2 catches the Green-identity sign, spheroid Jacobian, and incident phase.
rng(20260912, 'twister');
aspects = [2, 0.5];
volume_cases = cell(2, 2);
for case_index = 1:2
    aspect = aspects(case_index);
    a = aspect;
    b = 1;
    k = 8 / max(a, b);
    lambda = 2*pi/k;
    n_m = 1;
    n_p = 1.1;
    f = k^2 * (n_p^2 - 1);
    khat = [sin(pi/6), 0, cos(pi/6)];
    directions = randn(30, 3);
    directions = directions ./ vecnorm(directions, 2, 2);
    points = directions .* (1.2 + 3.8*rand(30, 1)) .* [b, b, a];
    % The planned 64/96 pair misses 1e-10 here; use the verified 96/144 pair.
    reference96 = born_volume_quadrature(a, b, k, f, khat, points, 96);
    reference144 = born_volume_quadrature(a, b, k, f, khat, points, 144);
    reference_change = max(abs(reference96 - reference144));
    fprintf('T2 a/b=%g: volume 96/144 change %.3e\n', aspect, reference_change);
    assert(reference_change <= 1e-10, 'T2 volume reference did not converge.');
    volume_cases(case_index,:) = {points, reference144};
end
for case_index = 1:2
    aspect = aspects(case_index);
    a = aspect;
    b = 1;
    lambda = 2*pi*max(a,b)/8;
    points = volume_cases{case_index,1};
    reference = volume_cases{case_index,2};
    uB = born_rytov_spheroid(lambda, n_p, n_m, a, b, pi/6, 0, ...
        points(:,1), points(:,2), points(:,3), struct('N_theta', 128));
    surface_error = max(abs(uB - reference));
    assert(surface_error <= 1e-8 * max(abs(reference)), ...
        'T2 surface Born disagrees with volume quadrature.');
    fprintf('T2 a/b=%g: surface relative error %.3e PASS\n', ...
        aspect, surface_error/max(abs(reference)));
end

% T3 catches far-point phase cancellation and the small-Q forward limit.
rng(20260913, 'twister');
for aspect = [2, 0.5]
    a = aspect;
    b = 1;
    k = 8 / max(a, b);
    lambda = 2*pi/k;
    f = k^2 * (1.1^2 - 1);
    khat = [sin(pi/6), 0, cos(pi/6)];
    in_plane = [-cos(pi/6), 0, sin(pi/6)];
    angles = (0:180)' * pi/180;
    directions = cos(angles)*khat + sin(angles)*in_plane;
    random_directions = randn(50, 3);
    random_directions = random_directions ./ vecnorm(random_directions, 2, 2);
    directions = [directions; random_directions];
    radius = 1e7 * max(a, b);
    points = radius * directions;
    reference = rgd_farfield(a, b, k, f, khat, directions);
    uB = born_rytov_spheroid(lambda, 1.1, 1, a, b, pi/6, 0, ...
        points(:,1), points(:,2), points(:,3), struct('N_theta', 64));
    amplitude = uB * radius * exp(-1i*k*radius);
    scale = abs(f) * (4*pi*a*b^2/3) / (4*pi);
    far_error = max(abs(amplitude - reference));
    assert(far_error <= 1e-5 * scale, 'T3 Rayleigh-Gans far field disagrees.');
    fprintf('T3 a/b=%g: scaled far-field error %.3e PASS\n', aspect, far_error/scale);
end

% T5 covers all size/aspect/incidence corners and start-order stability.
% Use the reviewed physical-clearance bound; elliptic clearance alone can be
% too small along an oblate/prolate object's short axis for default N_max.
rng(20260916, 'twister');
for x = [1, 30, 100]
    for aspect = [0.25, 1, 4]
        a = aspect;
        b = 1;
        k = x/max(a,b);
        for theta_inc = [0, pi/3]
            directions = randn(20,3);
            directions = directions ./ vecnorm(directions,2,2);
            clearance = 0.02*max(a,b)/min(a,b);
            radii = [0.2+(0.8-clearance)*rand(10,1); ...
                1+clearance+(3-clearance)*rand(10,1)];
            assert(all(abs(radii-1)*min(a,b) >= 0.02*max(a,b)), ...
                'T5 point generation violated the physical-clearance bound.');
            points = directions.*radii.*[b,b,a];
            args = {2*pi/k, 1.01, 1, a, b, theta_inc, 0, points(:,1), points(:,2), points(:,3)};
            [uB, ~, ~, info] = born_rytov_spheroid(args{:});
            [uB_doubled, ~, ~, doubled_info] = born_rytov_spheroid(args{:}, ...
                struct('N_theta', 2*(ceil(2*k*max(a,b))+16)));
            rho = sqrt(sum(points.^2,2));
            scale = max(abs(uB_doubled), abs(info.f)*info.V./(4*pi*max(rho,max(a,b))));
            change = max(abs(uB-uB_doubled)./scale);
            assert(info.validated && doubled_info.validated, 'T5 envelope case did not converge.');
            assert(change <= 1e-8, 'T5 doubled starting order changed Born beyond tolerance.');
            assert(all(info.err(:) <= 1e-8) && ~any(info.unconverged(:)), ...
                'T5 converged metadata disagrees with validation.');
            fprintf('T5 x=%g, a/b=%g, incidence=%g: start-order change %.3e PASS\n', ...
                x, aspect, theta_inc*180/pi, change);
        end
    end
end
[near_uB, ~, ~, near_info] = born_rytov_spheroid(2*pi/30, 1.01, 1, 1, 1, 0, 0, 1.02, 0, 0);
assert(near_info.validated && near_info.N_theta <= 1200, 'T5 r/a=1.02 did not converge.');
[group_uB, ~, ~, group_info] = born_rytov_spheroid(2*pi/30, 1.01, 1, 1, 1, ...
    0, 0, [1.02,3], [0,0], [0,0]);
assert(group_info.validated && group_info.N_theta(1)>group_info.N_theta(2), ...
    'T5 points were not accepted at their individual convergence orders.');
assert(abs(group_uB(1)-near_uB)<=1e-12*abs(near_uB), 'T5 grouped refinement changed the near field.');
X = reshape(linspace(1.2,3,42),7,3,2);
[uB, uR, u0, info] = born_rytov_spheroid(2*pi, 1.01, 1, 1, 1, 0, 0, X, 0*X, 0*X);
[chunked_uB, ~, ~, chunked_info] = born_rytov_spheroid(2*pi, 1.01, 1, 1, 1, ...
    0, 0, X, 0*X, 0*X, struct('chunk',1));
assert(chunked_info.validated && max(abs(uB(:)-chunked_uB(:)))<=1e-12*max(abs(uB(:))), ...
    'T5 chunking changed the adaptive Born result.');
arrays = {uB,uR,u0,info.N_theta,info.err,info.inside,info.unconverged};
assert(all(cellfun(@(value) isequal(size(value),size(X)), arrays)), 'T5 output shape was lost.');
assert(info.validated, 'T5 shaped inputs did not validate.');
[scalar_uB, scalar_uR, scalar_u0, scalar_info] = born_rytov_spheroid(2*pi, 1.01, 1, 1, 1, 0, 0, 2, 0, 0);
assert(all(cellfun(@isscalar, {scalar_uB,scalar_uR,scalar_u0,scalar_info.N_theta, ...
    scalar_info.err,scalar_info.inside,scalar_info.unconverged})), 'T5 scalar shape was lost.');
fprintf('T5 near-surface and shape checks PASS\n');
fprintf('test_born_rytov: T1/T2/T3/T4/T5 PASS (%.2f s)\n', toc(test_started));

function caught = assert_error(callback, identifier)
try
    callback();
catch caught
    assert(strcmp(caught.identifier,identifier), 'Expected %s, got %s.', identifier, caught.identifier);
    return
end
error('test_born_rytov:ExpectedError', 'Expected %s; failure path not exercised.', identifier);
end

function [u, extended] = born_sphere_series(a, k, f, khat, points, lmax, dephased)
% Addition-theorem reference: Bessel radial functions and independent integrals.
rho = sqrt(sum(points.^2,2));
mu = zeros(size(rho));
mu(rho>0) = points(rho>0,:)*khat.'./rho(rho>0);
extended = complex(zeros(size(rho)));
previous = zeros(size(rho));
legendre = ones(size(rho));
for ell = 0:lmax+10
    full_radial = integral(@(r) r.^2.*spherical_j(ell,k*r).^2, 0, a, ...
        'RelTol',1e-12,'AbsTol',0);
    radial = complex(zeros(size(rho)));
    outside = rho>=a;
    if dephased
        assert(all(outside) && all(k*rho>10*(lmax+10)^2), 'Far Hankel sum used outside its stable range.');
        z = k*rho;
        term = 1i^(-ell-1)./z;
        h = term;
        for m = 1:ell
            term = term.*(1i*(ell+m)*(ell-m+1)/(2*m))./z;
            h = h+term;
        end
        radial = h*full_radial;
    else
        radial(outside) = spherical_h(ell,k*rho(outside))*full_radial;
        for j = find(rho>0 & ~outside).'
            lower = integral(@(r) r.^2.*spherical_j(ell,k*r).^2, 0, rho(j), ...
                'RelTol',1e-12,'AbsTol',0);
            upper = integral(@(r) radial_jh(r,ell,k), rho(j), a, 'RelTol',1e-12,'AbsTol',0);
            radial(j) = spherical_h(ell,k*rho(j))*lower+spherical_j(ell,k*rho(j))*upper;
        end
        if ell == 0
            radial(rho==0) = integral(@(r) radial_jh(r,0,k), 0, a, 'RelTol',1e-12,'AbsTol',0);
        end
    end
    extended = extended+f*1i*k*(2*ell+1)*1i^ell*legendre.*radial;
    if ell == lmax
        u = extended;
    end
    following = ((2*ell+1)*mu.*legendre-ell*previous)/(ell+1);
    previous = legendre;
    legendre = following;
end
end

function value = spherical_j(ell,z)
value = sqrt(pi./(2*z)).*besselj(ell+0.5,z);
value(z==0) = double(ell==0);
end

function value = spherical_h(ell,z)
value = sqrt(pi./(2*z)).*(besselj(ell+0.5,z)+1i*bessely(ell+0.5,z));
end

function value = radial_jh(r,ell,k)
value = r.^2.*spherical_j(ell,k*r).*spherical_h(ell,k*r);
value(r==0) = 0;
end

function u = naive_born_surface(a, k, f, points, N)
% Test-only sphere surface rule with the intentionally unstable exp(ikR).
[g,wg] = gauss_legendre(N);
theta = pi*(g+1)/2;
phi = (0:2*N-1)*(pi/N);
[T,P] = ndgrid(theta,phi);
weights = (pi*wg/2)*ones(1,2*N)*(pi/N);
xyz = a*[sin(T(:)).*cos(P(:)),sin(T(:)).*sin(P(:)),cos(T(:))];
normal = xyz.*(a*sin(T(:)).*weights(:));
incident = exp(1i*k*xyz(:,3));
w = xyz(:,3).*incident/(2i*k);
dnw = normal(:,3).*incident.*(1+1i*k*xyz(:,3))/(2i*k);
u = complex(zeros(size(points,1),1));
for j = 1:size(points,1)
    offset = xyz-points(j,:);
    R = sqrt(sum(offset.^2,2));
    G = exp(1i*k*R)./(4*pi*R);
    dnG = (1i*k-1./R).*G.*sum(offset.*normal,2)./R;
    u(j) = f*(G.'*dnw-dnG.'*w);
end
end

function u = born_volume_quadrature(a, b, k, f, khat, points, N)
% Independent volume integral: GL radius and polar angle, trapezoid azimuth.
[g, wg] = gauss_legendre(N);
s = (g+1)/2;
ws = wg/2;
theta = pi*(g+1)/2;
wt = pi*wg/2;
phi = (0:N-1)' * (2*pi/N);
[S, T, P] = ndgrid(s, theta, phi);
[WS, WT, ~] = ndgrid(ws, wt, phi);
sources = [b*S(:).*sin(T(:)).*cos(P(:)), ...
           b*S(:).*sin(T(:)).*sin(P(:)), a*S(:).*cos(T(:))];
weights = a*b^2 * S(:).^2 .* sin(T(:)) .* WS(:) .* WT(:) * (2*pi/N);
density = weights .* exp(1i*k*(sources*khat.'));
u = complex(zeros(size(points,1), 1));
for j = 1:size(points,1)
    R = sqrt(sum((points(j,:) - sources).^2, 2));
    u(j) = f * sum(exp(1i*k*R) ./ (4*pi*R) .* density);
end
end

function amplitude = rgd_farfield(a, b, k, f, khat, directions)
q = k*(directions - khat);
Q = sqrt(b^2*(q(:,1).^2 + q(:,2).^2) + a^2*q(:,3).^2);
form = ones(size(Q));
small = Q < 1e-3;
form(small) = 1 - Q(small).^2/10;
form(~small) = 3*(sin(Q(~small)) - Q(~small).*cos(Q(~small))) ./ Q(~small).^3;
amplitude = f * (4*pi*a*b^2/3) / (4*pi) * form;
end
