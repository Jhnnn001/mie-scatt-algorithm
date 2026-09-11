function eg = sph_eigen(type, mu, c, L)
%SPH_EIGEN Angular spheroidal eigenvalues and Legendre coefficients.

narginchk(4, 4);
type = checked_type(type);
if ~(isnumeric(mu) && isreal(mu) && isscalar(mu) && isfinite(mu) && ...
        mu >= 0 && mu == fix(mu) && mu <= 300)
    error('sph_eigen:InvalidMu', 'mu must be an integer in [0, 300].');
end
if ~(isnumeric(c) && isscalar(c) && isfinite(real(c)) && isfinite(imag(c)))
    error('sph_eigen:InvalidC', 'c must be a finite numeric scalar.');
end
if ~(isnumeric(L) && isreal(L) && isscalar(L) && isfinite(L) && ...
        L >= mu && L == fix(L) && L <= 400)
    error('sph_eigen:InvalidL', 'L must be an integer in [mu, 400].');
end
mu = double(mu);
c = double(c);
L = double(L);

rmax = (L - mu) + ceil(abs(c)) + 40;
while true
    if rmax > 4000
        error('sph_eigen:NoConvergence', ...
            'Legendre coefficients require rmax above the limit of 4000.');
    end
    base = solve_truncation(type, mu, c, L, rmax);
    next_rmax = max(rmax + 2, ceil(1.5 * rmax));
    if next_rmax > 4000
        error('sph_eigen:NoConvergence', ...
            'Legendre coefficients require rmax above the limit of 4000.');
    end
    larger = solve_truncation(type, mu, c, L, next_rmax);
    [d, lambda, metrics] = align_retained(base, larger, mu, L);
    condition = cond(d);
    if ~isfinite(condition) || condition >= 1e8
        error('sph_eigen:ExceptionalPoint', ...
            'Retained eigenvector condition number %.3g is not below 1e8.', ...
            condition);
    end
    strict_ok = metrics.difference < 1e-13;
    floor_ok = metrics.difference <= 3e-13 && ...
        metrics.common_difference <= 3e-13 && metrics.tail <= 1e-14 && ...
        metrics.lambda_difference <= 100 * eps && ...
        metrics.backward_error <= 100 * eps;
    if strict_ok || floor_ok
        eg = struct('type', type, 'mu', mu, 'c', c, 'L', L, ...
            'l', (mu:L).', 'lambda', lambda, 'd', d, ...
            'r', (0:next_rmax).', 'rmax', next_rmax, ...
            'cond', condition, ...
            'rmax_delta', metrics.difference, ...
            'rmax_tail', metrics.tail, ...
            'rmax_floor_used', ~strict_ok && floor_ok);
        return
    end
    if metrics.tail <= 1e-14 && metrics.backward_error <= 100 * eps && ...
            (metrics.difference > 3e-13 || ...
            metrics.lambda_difference > 100 * eps)
        error('sph_eigen:EigenSensitivity', ...
            ['Backward-stable eigenpairs remain truncation-sensitive ' ...
            '(delta d %.3g, relative delta lambda %.3g, cond %.3g).'], ...
            metrics.difference, metrics.lambda_difference, condition);
    end
    rmax = next_rmax;
end
end

function result = solve_truncation(type, mu, c, L, rmax)
r = (0:rmax).';
c2 = c^2;
if strcmp(type, 'oblate')
    c2 = -c2;
end

blocks = repmat(struct('rows', [], 'lambda', [], 'vectors', [], ...
    'backward', []), 2, 1);
for parity = 0:1
    rows = find(mod(r, 2) == parity);
    rp = r(rows);
    ell = mu + rp;
    diagonal = ell .* (ell + 1) + ...
        (2 * ell .* (ell + 1) - 2 * mu^2 - 1) * c2 ./ ...
        ((2 * ell - 1) .* (2 * ell + 3));
    coupling_r = rp(1:end-1);
    off_diagonal = c2 * sqrt((2 * mu + coupling_r + 2) .* ...
        (2 * mu + coupling_r + 1) .* (coupling_r + 2) .* ...
        (coupling_r + 1)) ./ ((2 * mu + 2 * coupling_r + 3) .* ...
        sqrt((2 * mu + 2 * coupling_r + 1) .* ...
        (2 * mu + 2 * coupling_r + 5)));
    matrix = diag(diagonal) + diag(off_diagonal, 1) + ...
        diag(off_diagonal, -1);
    [vectors, lambda] = eig(matrix, 'vector');
    [~, order] = sortrows([real(lambda), imag(lambda)], [1, 2]);
    lambda = lambda(order);
    vectors = vectors(:, order);
    backward = zeros(size(lambda));
    matrix_norm = norm(matrix, inf);
    for k = 1:size(vectors, 2)
        bilinear_norm = vectors(:, k).' * vectors(:, k);
        if abs(bilinear_norm) <= eps * norm(vectors(:, k))^2
            error('sph_eigen:ExceptionalPoint', ...
                'An eigenvector has a vanishing bilinear norm.');
        end
        vectors(:, k) = vectors(:, k) / sqrt(bilinear_norm);
        target = min(k, size(vectors, 1));
        value = vectors(target, k);
        if real(value) < 0 || (real(value) == 0 && imag(value) < 0)
            vectors(:, k) = -vectors(:, k);
        end
        residual = matrix * vectors(:, k) - lambda(k) * vectors(:, k);
        scale = (matrix_norm + abs(lambda(k))) * norm(vectors(:, k));
        backward(k) = norm(residual) / max(realmin, scale);
    end
    blocks(parity + 1) = struct('rows', rows, 'lambda', lambda, ...
        'vectors', vectors, 'backward', backward);
end

d = retained_vectors(blocks, rmax, mu, L);
condition = cond(d);
if ~isfinite(condition) || condition >= 1e8
    error('sph_eigen:ExceptionalPoint', ...
        'Retained eigenvector condition number %.3g is not below 1e8.', ...
        condition);
end
result = struct('rmax', rmax, 'blocks', blocks, ...
    'condition', condition, 'd', d, ...
    'lambda', retained_eigenvalues(blocks, mu, L));
end

function d = retained_vectors(blocks, rmax, mu, L)
d = zeros(rmax + 1, L - mu + 1);
for parity = 0:1
    columns = find(mod((mu:L) - mu, 2) == parity);
    count = numel(columns);
    if count > 0
        block = blocks(parity + 1);
        d(block.rows, columns) = block.vectors(:, 1:count);
    end
end
end

function lambda = retained_eigenvalues(blocks, mu, L)
lambda = zeros(L - mu + 1, 1, 'like', blocks(1).lambda);
for parity = 0:1
    columns = find(mod((mu:L) - mu, 2) == parity);
    count = numel(columns);
    if count > 0
        lambda(columns) = blocks(parity + 1).lambda(1:count);
    end
end
end

function [d, lambda, metrics] = align_retained(base, larger, mu, L)
d = zeros(larger.rmax + 1, L - mu + 1);
lambda = zeros(L - mu + 1, 1, 'like', larger.blocks(1).lambda);
difference = 0;
common_difference = 0;
tail = 0;
lambda_difference = 0;
backward_error = 0;
for parity = 0:1
    columns = find(mod((mu:L) - mu, 2) == parity);
    if isempty(columns)
        continue
    end
    base_rows = base.blocks(parity + 1).rows;
    larger_block = larger.blocks(parity + 1);
    common_count = numel(base_rows);
    base_vectors = base.d(base_rows, columns);
    candidate_common = larger_block.vectors(1:common_count, :);
    base_norms = sqrt(sum(abs(base_vectors).^2, 1)).';
    candidate_norms = sqrt(sum(abs(candidate_common).^2, 1));
    overlap = abs(base_vectors' * candidate_common) ./ ...
        (base_norms * candidate_norms);
    base_lambda = base.lambda(columns);
    candidate_lambda = larger_block.lambda;
    eigenvalue_scale = max(1, max(abs(base_lambda), abs(candidate_lambda).'));
    eigenvalue_distance = abs(base_lambda - candidate_lambda.') ./ ...
        eigenvalue_scale;
    matches = greedy_match(eigenvalue_distance, overlap);
    for j = 1:numel(columns)
        column = columns(j);
        vector = zeros(larger.rmax + 1, 1);
        vector(larger_block.rows) = larger_block.vectors(:, matches(j));
        padded_base = zeros(larger.rmax + 1, 1);
        padded_base(1:base.rmax + 1) = base.d(:, column);
        if norm(vector + padded_base) < norm(vector - padded_base)
            vector = -vector;
        end
        difference = max(difference, norm(vector - padded_base));
        common_difference = max(common_difference, ...
            norm(vector(1:base.rmax + 1) - base.d(:, column)));
        tail = max(tail, norm(vector(base.rmax + 2:end)));
        target = column;
        value = vector(target);
        if real(value) < 0 || (real(value) == 0 && imag(value) < 0)
            vector = -vector;
        end
        d(:, column) = vector;
        lambda(column) = larger_block.lambda(matches(j));
        lambda_scale = max([1, abs(lambda(column)), abs(base.lambda(column))]);
        lambda_difference = max(lambda_difference, ...
            abs(lambda(column) - base.lambda(column)) / lambda_scale);
        backward_error = max([backward_error, ...
            base.blocks(parity + 1).backward(j), ...
            larger_block.backward(matches(j))]);
    end
end
metrics = struct('difference', difference, ...
    'common_difference', common_difference, 'tail', tail, ...
    'lambda_difference', lambda_difference, ...
    'backward_error', backward_error);
end

function matches = greedy_match(distance, overlap)
count = size(distance, 1);
matches = zeros(count, 1);
available_rows = true(count, 1);
available_columns = true(size(distance, 2), 1);
for k = 1:count
    remaining_distance = distance;
    remaining_distance(~available_rows, :) = Inf;
    remaining_distance(:, ~available_columns) = Inf;
    minimum = min(remaining_distance(:));
    tied = remaining_distance <= minimum + 100 * eps;
    scores = overlap;
    scores(~tied) = -Inf;
    [~, linear_index] = max(scores(:));
    [row, column] = ind2sub(size(scores), linear_index);
    matches(row) = column;
    available_rows(row) = false;
    available_columns(column) = false;
end
end

function type = checked_type(type)
if isstring(type) && isscalar(type)
    type = char(type);
end
if ~(ischar(type) && isrow(type))
    error('sph_eigen:InvalidType', ...
        'type must be ''prolate'' or ''oblate''.');
end
type = lower(type);
if ~strcmp(type, 'prolate') && ~strcmp(type, 'oblate')
    error('sph_eigen:InvalidType', ...
        'type must be ''prolate'' or ''oblate''.');
end
end
