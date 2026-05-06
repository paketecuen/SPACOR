function result = sp1_identify_bh_loop_parameters(data, opts)
%SP1_IDENTIFY_BH_LOOP_PARAMETERS Identify R and B-H loop parameters.
%
% For each candidate R:
%
%   lambda_R(t) = integral(v - R i) dt
%
% and the B-H relation is fitted as:
%
%   i = a1 lambda + a3 lambda^3 + h b(lambda, dlambda/dt).
%
% Returned physical parameters:
%
%   R
%   L0 = 1/a1
%   saturationCubic = a3
%   hysteresisCurrent = h

if nargin < 2
    opts = struct();
end

[workData, preconditioning] = sp1_precondition_magnetic_measurements(data, opts);
t = workData.t(:);
i = workData.i(:);
v = workData.v(:);
rGrid = nv3_get_option(opts, 'rGrid', linspace(0, 0.6, 301));
includeHysteresis = nv3_get_option(opts, 'includeHysteresis', true);
hysteresisSmoothing = nv3_get_option(opts, 'hysteresisSmoothing', ...
    default_hysteresis_smoothing(data));
closureWeight = nv3_get_option(opts, 'closureWeight', 0.02);
passivityWeight = nv3_get_option(opts, 'passivityWeight', 10);

rows = cell(numel(rGrid), 13);
best = struct('score', Inf);
for idx = 1:numel(rGrid)
    R = rGrid(idx);
    lambda = reconstruct_lambda(t, v, i, R);
    dlambda = v - R * i;
    [X, basis] = bh_basis(lambda, dlambda, includeHysteresis, ...
        hysteresisSmoothing);
    [coef, fit] = normalized_solve(X, i);
    iHat = X * coef;
    residual = norm(i - iHat) / max(norm(i), eps);
    closure = abs(lambda(end) - lambda(1)) / ...
        max(max(lambda) - min(lambda), eps);
    [L0, satCubic, hCurrent] = physical_parameters(coef, includeHysteresis);
    passivePenalty = passivity_penalty(R, L0, hCurrent);
    score = residual + closureWeight * closure + passivityWeight * passivePenalty;
    magneticLoss = magnetic_loss_per_cycle(t, i, dlambda, data.f0);

    rows(idx, :) = {R, score, residual, closure, fit.conditionNumber, ...
        fit.rank, coef(:).', L0, satCubic, hCurrent, magneticLoss, ...
        passivePenalty, basis.lambdaScale};

    if score < best.score
        best = struct('R', R, 'score', score, 'residual', residual, ...
            'closure', closure, 'conditionNumber', fit.conditionNumber, ...
            'rank', fit.rank, 'coef', coef, 'lambda', lambda, ...
            'dlambda', dlambda, 'iHat', iHat, 'basis', basis, ...
            'L0', L0, 'saturationCubic', satCubic, ...
            'hysteresisCurrent', hCurrent, ...
            'magneticLossPerCycle', magneticLoss, ...
            'passivePenalty', passivePenalty);
    end
end

ranking = cell2table(rows, 'VariableNames', { ...
    'R', 'Score', 'Residual', 'Closure', 'ConditionNumber', 'Rank', ...
    'Coefficients', 'L0', 'SaturationCubic', 'HysteresisCurrent', ...
    'MagneticLossPerCycle', 'PassivePenalty', 'LambdaScale'});

truth = truth_comparison(data, best);

result = struct();
result.model = 'bh_loop_parametric';
result.ranking = ranking;
result.R = best.R;
result.L0 = best.L0;
result.saturationCubic = best.saturationCubic;
result.hysteresisCurrent = best.hysteresisCurrent;
result.score = best.score;
result.residual = best.residual;
result.closure = best.closure;
result.conditionNumber = best.conditionNumber;
result.rank = best.rank;
result.coefficients = best.coef;
result.lambda = best.lambda;
result.dlambda = best.dlambda;
result.iFit = best.iHat;
result.basis = best.basis;
result.magneticLossPerCycle = best.magneticLossPerCycle;
result.truth = truth;
result.preconditioning = preconditioning;
result.options = struct('rGrid', rGrid, ...
    'includeHysteresis', includeHysteresis, ...
    'hysteresisSmoothing', hysteresisSmoothing);

end

function lambda = reconstruct_lambda(t, v, i, R)
lambda = cumtrapz(t, v - R * i);
lambda = lambda - mean(lambda, 'omitnan');
end

function [X, basis] = bh_basis(lambda, dlambda, includeHysteresis, smoothing)
lambdaScale = max(abs(lambda));
lambdaScale = max(lambdaScale, eps);
u = lambda ./ lambdaScale;
branchScale = max(abs(dlambda));
branchScale = max(smoothing * branchScale, eps);
branch = tanh(dlambda ./ branchScale);
shape = 0.25 + 0.75 * sqrt(abs(u));
X = [lambda, lambda .^ 3];
if includeHysteresis
    X = [X, shape .* branch];
end
basis = struct('lambdaScale', lambdaScale, 'branch', branch, ...
    'shape', shape, 'hysteresisBasis', shape .* branch);
end

function [coef, fit] = normalized_solve(X, y)
mask = all(isfinite(X), 2) & isfinite(y);
Xv = X(mask, :);
yv = y(mask);
m = size(X, 2);
coef = NaN(m, 1);
fit = struct('rank', 0, 'conditionNumber', Inf, ...
    'singularValues', NaN);
if size(Xv, 1) < m
    return;
end
colScale = sqrt(mean(Xv .^ 2, 1)).';
yScale = sqrt(mean(yv .^ 2));
colScale(~isfinite(colScale) | colScale <= eps) = 1;
yScale = max(yScale, eps);
Xs = Xv ./ colScale.';
ys = yv ./ yScale;
sv = svd(Xs, 'econ');
fit.singularValues = sv(:).';
if isempty(sv) || max(sv) <= eps
    return;
end
tol = max(size(Xs)) * eps(max(sv));
fit.rank = sum(sv > tol);
if numel(sv) >= m && sv(m) > tol
    fit.conditionNumber = sv(1) / sv(m);
end
if fit.rank < m
    return;
end
gamma = Xs \ ys;
coef = gamma .* yScale ./ colScale;
end

function [L0, satCubic, hCurrent] = physical_parameters(coef, includeHysteresis)
L0 = 1 / coef(1);
satCubic = coef(2);
if includeHysteresis && numel(coef) >= 3
    hCurrent = coef(3);
else
    hCurrent = 0;
end
end

function penalty = passivity_penalty(R, L0, hCurrent)
penalty = 0;
if ~isfinite(R) || R < 0
    penalty = penalty + 1;
end
if ~isfinite(L0) || L0 <= 0
    penalty = penalty + 1;
end
if ~isfinite(hCurrent) || hCurrent < 0
    penalty = penalty + 0.1;
end
end

function loss = magnetic_loss_per_cycle(t, i, dlambda, f0)
samplesPerCycle = max(4, round((1 / f0) / median(diff(t))));
nCycles = floor(numel(t) / samplesPerCycle);
losses = zeros(nCycles, 1);
for idx = 1:nCycles
    range = (idx - 1) * samplesPerCycle + (1:samplesPerCycle);
    losses(idx) = trapz(t(range), i(range) .* dlambda(range));
end
loss = median(losses, 'omitnan');
end

function smoothing = default_hysteresis_smoothing(data)
if isfield(data, 'truth') && isfield(data.truth, 'hysteresisSmoothing')
    smoothing = data.truth.hysteresisSmoothing;
else
    smoothing = 0.04;
end
end

function truth = truth_comparison(data, best)
truth = struct();
if ~isfield(data, 'truth')
    return;
end
fields = {'R', 'L0', 'saturationCubic', 'hysteresisCurrent'};
est = [best.R, best.L0, best.saturationCubic, best.hysteresisCurrent];
for idx = 1:numel(fields)
    field = fields{idx};
    if isfield(data.truth, field)
        truth.([field 'Truth']) = data.truth.(field);
        truth.([field 'RelativeError']) = abs(est(idx) - data.truth.(field)) / ...
            max(abs(data.truth.(field)), eps);
    end
end
if isfield(data, 'clean') && isfield(data.clean, 'lambda')
    lambdaTruth = data.clean.lambda(:) - mean(data.clean.lambda(:), 'omitnan');
    lambdaHat = best.lambda(:);
    scale = dot(lambdaHat, lambdaTruth) / max(dot(lambdaHat, lambdaHat), eps);
    truth.LambdaRelativeError = norm(scale * lambdaHat - lambdaTruth) / ...
        max(norm(lambdaTruth), eps);
end
if isfield(data, 'energy') && isfield(data.energy, 'MedianMagneticLossPerCycle')
    truth.MagneticLossTruth = data.energy.MedianMagneticLossPerCycle;
    truth.MagneticLossRelativeError = abs(best.magneticLossPerCycle - ...
        data.energy.MedianMagneticLossPerCycle) / ...
        max(abs(data.energy.MedianMagneticLossPerCycle), eps);
end
end
