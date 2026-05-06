function result = sp1_identify_flux_curve(data, ~, opts)
%SP1_IDENTIFY_FLUX_CURVE Estimate lambda=f(i) from v-Ri.
%
% The resistance is selected by minimizing the single-valuedness error of
% lambda(i) over an odd polynomial basis.

if nargin < 3
    opts = struct();
end

degree = nv3_get_option(opts, 'degree', 9);
rGrid = nv3_get_option(opts, 'rGrid', linspace(0.2, 1.2, 241));
if mod(degree, 2) == 0
    degree = degree + 1;
end
orders = 1:2:degree;

[workData, preconditioning] = sp1_precondition_magnetic_measurements(data, opts);
t = workData.t(:);
i = workData.i(:);
v = workData.v(:);
rows = cell(numel(rGrid), 5);
best = struct('score', Inf);
for idx = 1:numel(rGrid)
    R = rGrid(idx);
    lambda = cumtrapz(t, v - R * i);
    lambda = lambda - mean(lambda);
    [coef, fitLambda, score, basisCondition] = fit_odd_lambda(i, lambda, orders);
    rows(idx, :) = {R, score, basisCondition, norm(lambda), mat2str(coef(:).')};
    if score < best.score
        best = struct('R', R, 'score', score, 'basisCondition', basisCondition, ...
            'lambdaRaw', lambda, 'lambdaFit', fitLambda, 'coef', coef);
    end
end

lambdaError = truth_lambda_error(data, best.lambdaFit);

currentGrid = linspace(min(i), max(i), 500).';
B = odd_basis(currentGrid, orders);
lambdaGrid = B * best.coef;
LincGrid = odd_basis_derivative(currentGrid, orders) * best.coef;

result = struct();
result.model = 'flux_curve_constitutive';
result.ranking = cell2table(rows, 'VariableNames', ...
    {'R', 'SingleValuednessError', 'BasisCondition', 'LambdaNorm', 'Coefficients'});
result.R = best.R;
result.score = best.score;
result.basisCondition = best.basisCondition;
result.lambdaRaw = best.lambdaRaw;
result.lambdaFit = best.lambdaFit;
result.lambdaRelativeError = lambdaError;
result.coefficients = best.coef;
result.orders = orders;
result.preconditioning = preconditioning;
[lambdaTruthGrid, lincTruthGrid] = truth_curve(data, currentGrid);
result.curve = table(currentGrid, lambdaGrid, LincGrid, ...
    lambdaTruthGrid, lincTruthGrid, ...
    'VariableNames', {'Current', 'LambdaEstimate', 'LincEstimate', ...
    'LambdaTruth', 'LincTruth'});

end

function [lambda, linc] = truth_curve(data, i)
if isfield(data, 'truth') && isfield(data.truth, 'law') && ...
        strcmpi(char(data.truth.law), 'current_saturation')
    u = i ./ data.truth.Is;
    L = data.truth.Lmin + ...
        (data.truth.Lsat0 - data.truth.Lmin) ./ cosh(u) .^ 2;
    dLdi = -2 * (data.truth.Lsat0 - data.truth.Lmin) .* ...
        tanh(u) ./ (data.truth.Is .* cosh(u) .^ 2);
    lambda = L .* i;
    linc = L + i .* dLdi;
elseif isfield(data, 'truth') && isfield(data.truth, 'Lmin') && ...
        isfield(data.truth, 'L0') && isfield(data.truth, 'Is')
    lambda = sp1_flux_curve(i, data.truth);
    linc = sp1_incremental_inductance(i, data.truth);
else
    lambda = NaN(size(i));
    linc = NaN(size(i));
end
end

function lambdaError = truth_lambda_error(data, lambdaFit)
lambdaError = NaN;
if isfield(data, 'clean') && isfield(data.clean, 'lambda')
    truthLambda = data.clean.lambda(:);
    truthLambda = truthLambda - mean(truthLambda);
    offset = median(lambdaFit - truthLambda, 'omitnan');
    lambdaError = norm((lambdaFit - offset) - truthLambda) / ...
        max(norm(truthLambda), eps);
end
end

function [coef, fitLambda, score, basisCondition] = fit_odd_lambda(i, lambda, orders)
B = odd_basis(i, orders);
sv = svd(B, 'econ');
if isempty(sv) || min(sv) <= eps(max(sv))
    basisCondition = Inf;
else
    basisCondition = max(sv) / min(sv);
end
coef = B \ lambda;
fitLambda = B * coef;
score = norm(lambda - fitLambda) / max(norm(lambda), eps);
end

function B = odd_basis(i, orders)
i = i(:);
scale = max(abs(i));
scale = max(scale, eps);
u = i ./ scale;
B = zeros(numel(i), numel(orders));
for idx = 1:numel(orders)
    B(:, idx) = scale * u .^ orders(idx);
end
end

function B = odd_basis_derivative(i, orders)
i = i(:);
scale = max(abs(i));
scale = max(scale, eps);
u = i ./ scale;
B = zeros(numel(i), numel(orders));
for idx = 1:numel(orders)
    n = orders(idx);
    B(:, idx) = n * u .^ (n - 1);
end
end
