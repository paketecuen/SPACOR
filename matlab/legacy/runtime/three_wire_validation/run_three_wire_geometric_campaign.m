function campaign = run_three_wire_geometric_campaign(mode)
%RUN_THREE_WIRE_GEOMETRIC_CAMPAIGN Time-domain identifiability tests for 3-wire RLC.
%
% The campaign asks a strictly time-domain question:
%   Does the measured trajectory span enough independent geometric directions
%   to estimate the selected physical parameterization?
%
% The tested physical model is delta parallel RLC:
%   i_xy = G_xy v_xy + Gamma_xy integral(v_xy) + C_xy dv_xy/dt.
%
% Symmetry classes reduce the number of free unknowns:
%   balanced        : 3 free parameters
%   two_equal       : 6 free parameters, ab=bc and ca free
%   full_unbalanced : 9 free parameters

if nargin < 1 || isempty(mode)
    mode = 'full';
end

config = nv3_default_config();
if strcmpi(mode, 'quick')
    config.duration = 0.4;
    windowCyclesList = [0.5, 1.0];
else
    config.duration = 0.8;
    windowCyclesList = [0.25, 0.5, 1.0, 2.0];
end
config.fs = 20000;
config.f0 = 50;

outDir = config.resultsDir;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

t = (0:round(config.duration * config.fs) - 1).' / config.fs;
excitationNames = {'planar_balanced', 'planar_unbalanced', ...
    'weak_curvature', 'rich_curve', 'localized_transient'};
symmetryNames = {'balanced', 'two_equal', 'full_unbalanced'};

summaryRows = {};
parameterRows = {};
caseIndex = 0;

for e = 1:numel(excitationNames)
    excitationName = excitationNames{e};
    curve = build_excitation(excitationName, t, config.f0);
    for s = 1:numel(symmetryNames)
        symmetryName = symmetryNames{s};
        truth = build_truth(symmetryName);
        data = synthesize_delta_rlc(curve, truth);
        for wc = 1:numel(windowCyclesList)
            caseIndex = caseIndex + 1;
            windowCycles = windowCyclesList(wc);
            result = evaluate_case(data, truth, symmetryName, windowCycles, config);
            summaryRows(end + 1, :) = { ...
                caseIndex, excitationName, symmetryName, result.nParameters, ...
                windowCycles, result.coverage, result.medianRank, ...
                result.medianCondition, result.p90Condition, ...
                result.maxRelativeError, result.medianPowerResidual, ...
                result.pass}; %#ok<AGROW>
            parameterRows = [parameterRows; parameter_rows(caseIndex, ...
                excitationName, symmetryName, windowCycles, result)]; %#ok<AGROW>
        end
    end
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'CaseIndex', 'Excitation', 'Symmetry', 'NParameters', ...
    'WindowCycles', 'Coverage', 'MedianRank', 'MedianCondition', ...
    'P90Condition', 'MaxRelativeError', 'MedianPowerResidual', 'Pass'});

if isempty(parameterRows)
    parameters = table();
else
    parameters = cell2table(parameterRows, 'VariableNames', { ...
        'CaseIndex', 'Excitation', 'Symmetry', 'WindowCycles', ...
        'Parameter', 'Truth', 'MedianEstimate', 'RelativeError', 'IQR'});
end

prefix = 'three_wire_geometric_campaign';
if strcmpi(mode, 'quick')
    prefix = 'three_wire_geometric_quick';
end
writetable(summary, fullfile(outDir, [prefix '_summary.csv']));
writetable(parameters, fullfile(outDir, [prefix '_parameters.csv']));
save(fullfile(outDir, [prefix '.mat']), 'summary', 'parameters');
plot_summary(summary, outDir, prefix);

campaign = struct();
campaign.summary = summary;
campaign.parameters = parameters;
campaign.config = config;
campaign.outputPrefix = prefix;

disp(summary);
end

function curve = build_excitation(name, t, f0)
switch name
    case 'planar_balanced'
        vabTerms = struct('kind', {'cos'}, 'amp', {120}, 'mult', {1}, 'phase', {0});
        vbcTerms = struct('kind', {'cos'}, 'amp', {120}, 'mult', {1}, 'phase', {-2*pi/3});
    case 'planar_unbalanced'
        vabTerms = struct('kind', {'cos'}, 'amp', {120}, 'mult', {1}, 'phase', {0.22});
        vbcTerms = struct('kind', {'cos'}, 'amp', {83}, 'mult', {1}, 'phase', {-1.86});
    case 'weak_curvature'
        vabTerms = [ ...
            struct('kind', 'cos', 'amp', 120, 'mult', 1, 'phase', 0.00), ...
            struct('kind', 'cos', 'amp', 0.25, 'mult', 5, 'phase', 0.40)];
        vbcTerms = [ ...
            struct('kind', 'cos', 'amp', 95, 'mult', 1, 'phase', -2.05), ...
            struct('kind', 'cos', 'amp', 0.20, 'mult', 5, 'phase', -0.80)];
    case 'rich_curve'
        vabTerms = [ ...
            struct('kind', 'cos', 'amp', 120, 'mult', 1, 'phase', 0.00), ...
            struct('kind', 'cos', 'amp', 10, 'mult', 5, 'phase', 0.40), ...
            struct('kind', 'cos', 'amp', 5, 'mult', 7, 'phase', -1.00)];
        vbcTerms = [ ...
            struct('kind', 'cos', 'amp', 92, 'mult', 1, 'phase', -2.18), ...
            struct('kind', 'cos', 'amp', 7, 'mult', 5, 'phase', -0.75), ...
            struct('kind', 'cos', 'amp', 4, 'mult', 7, 'phase', 1.15)];
    case 'localized_transient'
        vabTerms = [ ...
            struct('kind', 'cos', 'amp', 120, 'mult', 1, 'phase', 0.00), ...
            struct('kind', 'gaussian', 'amp', 18, 'mult', 0, 'phase', 0.19)];
        vbcTerms = [ ...
            struct('kind', 'cos', 'amp', 92, 'mult', 1, 'phase', -2.18), ...
            struct('kind', 'gaussian', 'amp', -14, 'mult', 0, 'phase', 0.23)];
    otherwise
        error('Unknown excitation "%s".', name);
end

[vab, dvab, qab] = eval_terms(vabTerms, t, f0);
[vbc, dvbc, qbc] = eval_terms(vbcTerms, t, f0);
curve = struct();
curve.vab = vab;
curve.vbc = vbc;
curve.vca = -vab - vbc;
curve.dvab = dvab;
curve.dvbc = dvbc;
curve.dvca = -dvab - dvbc;
curve.qab = qab;
curve.qbc = qbc;
curve.qca = -qab - qbc;
end

function [v, dv, q] = eval_terms(terms, t, f0)
w0 = 2*pi*f0;
v = zeros(size(t));
dv = zeros(size(t));
q = zeros(size(t));
for k = 1:numel(terms)
    term = terms(k);
    switch term.kind
        case 'cos'
            w = term.mult * w0;
            arg = w * t + term.phase;
            v = v + term.amp * cos(arg);
            dv = dv - term.amp * w * sin(arg);
            q = q + term.amp / w * sin(arg);
        case 'gaussian'
            center = term.phase;
            sigma = 0.010;
            tau = (t - center) / sigma;
            g = exp(-(tau .^ 2));
            v = v + term.amp * g;
            dv = dv + term.amp * g .* (-2 * (t - center) / sigma^2);
            q = q + term.amp * sigma * sqrt(pi) / 2 * erf(tau);
        otherwise
            error('Unknown time term "%s".', term.kind);
    end
end
end

function truth = build_truth(symmetryName)
switch symmetryName
    case 'balanced'
        truth.G = [1/20, 1/20, 1/20];
        truth.Gamma = [1/0.080, 1/0.080, 1/0.080];
        truth.C = [80e-6, 80e-6, 80e-6];
    case 'two_equal'
        truth.G = [1/20, 1/20, 1/15];
        truth.Gamma = [1/0.080, 1/0.080, 1/0.060];
        truth.C = [80e-6, 80e-6, 60e-6];
    case 'full_unbalanced'
        truth.G = [1/20, 1/30, 1/15];
        truth.Gamma = [1/0.080, 1/0.120, 1/0.060];
        truth.C = [80e-6, 120e-6, 60e-6];
    otherwise
        error('Unknown symmetry "%s".', symmetryName);
end
truth.namesFull = {'Gab','Gbc','Gca','Gammaab','Gammabc','Gammaca','Cab','Cbc','Cca'};
truth.thetaFull = [truth.G(:); truth.Gamma(:); truth.C(:)];
end

function data = synthesize_delta_rlc(curve, truth)
iab = truth.G(1) * curve.vab + truth.Gamma(1) * curve.qab + truth.C(1) * curve.dvab;
ibc = truth.G(2) * curve.vbc + truth.Gamma(2) * curve.qbc + truth.C(2) * curve.dvbc;
ica = truth.G(3) * curve.vca + truth.Gamma(3) * curve.qca + truth.C(3) * curve.dvca;
data = curve;
data.iab = iab;
data.ibc = ibc;
data.ica = ica;
data.ia = iab - ica;
data.ib = ibc - iab;
data.ic = ica - ibc;
data.power = data.vab .* data.ia + data.vbc .* (data.ia + data.ib);
end

function result = evaluate_case(data, truth, symmetryName, windowCycles, config)
[S, freeNames, freeTruth] = symmetry_matrix(symmetryName, truth);
nParameters = numel(freeTruth);
windowSamples = max(12, round(windowCycles * config.fs / config.f0));
hopSamples = max(1, round(0.1 * config.fs / config.f0));
n = numel(data.vab);
starts = 1:hopSamples:(n - windowSamples + 1);

rankValues = NaN(numel(starts), 1);
conditionValues = Inf(numel(starts), 1);
powerResiduals = Inf(numel(starts), 1);
thetaValues = NaN(numel(starts), nParameters);
valid = false(numel(starts), 1);

for w = 1:numel(starts)
    idx = starts(w):(starts(w) + windowSamples - 1);
    [Xfull, y] = delta_full_system(data, idx);
    X = Xfull * S;
    [theta, fit] = normalized_solve(X, y);
    rankValues(w) = fit.rank;
    conditionValues(w) = fit.conditionNumber;
    if fit.rank == nParameters && isfinite(fit.conditionNumber)
        thetaValues(w, :) = theta.';
        fullTheta = S * theta;
        pred = predict_delta(data, idx, fullTheta);
        powerResiduals(w) = norm(data.power(idx) - pred.power) / max(norm(data.power(idx)), eps);
        valid(w) = fit.conditionNumber < config.conditionLimit && ...
            all(theta > 0) && powerResiduals(w) < 1e-9;
    end
end

medTheta = NaN(1, nParameters);
relErr = NaN(1, nParameters);
iqrTheta = NaN(1, nParameters);
for p = 1:nParameters
    values = thetaValues(valid, p);
    values = values(isfinite(values));
    if ~isempty(values)
        medTheta(p) = median(values);
        relErr(p) = abs(medTheta(p) - freeTruth(p)) / max(abs(freeTruth(p)), eps);
        iqrTheta(p) = percentile(values, 75) - percentile(values, 25);
    end
end

coverage = mean(valid);
result = struct();
result.nParameters = nParameters;
result.parameterNames = freeNames;
result.truth = freeTruth;
result.medianTheta = medTheta(:);
result.relativeError = relErr(:);
result.iqrTheta = iqrTheta(:);
result.coverage = coverage;
result.medianRank = median(rankValues(isfinite(rankValues)));
result.medianCondition = median_finite(conditionValues(valid));
result.p90Condition = percentile(conditionValues(valid), 90);
result.maxRelativeError = max(relErr, [], 'omitnan');
result.medianPowerResidual = median_finite(powerResiduals(valid));
result.pass = coverage >= 0.80 && result.maxRelativeError < 1e-8 && ...
    result.medianCondition < config.conditionLimit && result.medianPowerResidual < 1e-9;
end

function [S, names, theta] = symmetry_matrix(symmetryName, truth)
switch symmetryName
    case 'balanced'
        S = [
            1 0 0
            1 0 0
            1 0 0
            0 1 0
            0 1 0
            0 1 0
            0 0 1
            0 0 1
            0 0 1];
        names = {'G', 'Gamma', 'C'};
        theta = [truth.G(1); truth.Gamma(1); truth.C(1)];
    case 'two_equal'
        S = [
            1 0 0 0 0 0
            1 0 0 0 0 0
            0 1 0 0 0 0
            0 0 1 0 0 0
            0 0 1 0 0 0
            0 0 0 1 0 0
            0 0 0 0 1 0
            0 0 0 0 1 0
            0 0 0 0 0 1];
        names = {'G_ab_bc', 'G_ca', 'Gamma_ab_bc', 'Gamma_ca', 'C_ab_bc', 'C_ca'};
        theta = [truth.G(1); truth.G(3); truth.Gamma(1); truth.Gamma(3); truth.C(1); truth.C(3)];
    case 'full_unbalanced'
        S = eye(9);
        names = truth.namesFull;
        theta = truth.thetaFull;
    otherwise
        error('Unknown symmetry "%s".', symmetryName);
end
end

function [X, y] = delta_full_system(data, idx)
vab = data.vab(idx);
vbc = data.vbc(idx);
vca = data.vca(idx);
qab = data.qab(idx);
qbc = data.qbc(idx);
qca = data.qca(idx);
dvab = data.dvab(idx);
dvbc = data.dvbc(idx);
dvca = data.dvca(idx);
z = zeros(size(vab));
X = [
    vab, z, -vca, qab, z, -qca, dvab, z, -dvca
    -vab, vbc, z, -qab, qbc, z, -dvab, dvbc, z
    ];
y = [data.ia(idx); data.ib(idx)];
end

function [theta, fit] = normalized_solve(X, y)
mask = all(isfinite(X), 2) & isfinite(y);
Xv = X(mask, :);
yv = y(mask);
m = size(X, 2);
theta = NaN(m, 1);
fit = struct('rank', 0, 'conditionNumber', Inf, 'relativeResidual', Inf);
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
tol = max(size(Xs)) * eps(max(sv));
fit.rank = sum(sv > tol);
if numel(sv) >= m && sv(m) > tol
    fit.conditionNumber = sv(1) / sv(m);
end
if fit.rank < m
    return;
end
gamma = Xs \ ys;
theta = gamma .* yScale ./ colScale;
fit.relativeResidual = norm(yv - Xv * theta) / max(norm(yv), eps);
end

function pred = predict_delta(data, idx, thetaFull)
Gab = thetaFull(1);
Gbc = thetaFull(2);
Gca = thetaFull(3);
Gammaab = thetaFull(4);
Gammabc = thetaFull(5);
Gammaca = thetaFull(6);
Cab = thetaFull(7);
Cbc = thetaFull(8);
Cca = thetaFull(9);
iab = Gab * data.vab(idx) + Gammaab * data.qab(idx) + Cab * data.dvab(idx);
ibc = Gbc * data.vbc(idx) + Gammabc * data.qbc(idx) + Cbc * data.dvbc(idx);
ica = Gca * data.vca(idx) + Gammaca * data.qca(idx) + Cca * data.dvca(idx);
ia = iab - ica;
ib = ibc - iab;
pred.ia = ia;
pred.ib = ib;
pred.power = data.vab(idx) .* ia + data.vbc(idx) .* (ia + ib);
end

function rows = parameter_rows(caseIndex, excitationName, symmetryName, windowCycles, result)
rows = cell(numel(result.parameterNames), 9);
for p = 1:numel(result.parameterNames)
    rows(p, :) = {caseIndex, excitationName, symmetryName, windowCycles, ...
        result.parameterNames{p}, result.truth(p), result.medianTheta(p), ...
        result.relativeError(p), result.iqrTheta(p)};
end
end

function plot_summary(summary, outDir, prefix)
symmetryNames = unique(summary.Symmetry, 'stable');
excitationNames = unique(summary.Excitation, 'stable');
for wc = unique(summary.WindowCycles).'
    maskW = summary.WindowCycles == wc;
    coverage = nan(numel(excitationNames), numel(symmetryNames));
    condition = nan(numel(excitationNames), numel(symmetryNames));
    for e = 1:numel(excitationNames)
        for s = 1:numel(symmetryNames)
            mask = maskW & strcmp(summary.Excitation, excitationNames{e}) & ...
                strcmp(summary.Symmetry, symmetryNames{s});
            if any(mask)
                coverage(e, s) = summary.Coverage(mask);
                condition(e, s) = summary.MedianCondition(mask);
            end
        end
    end
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 450]);
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    imagesc(coverage, [0 1]);
    colorbar;
    title(sprintf('Valid coverage | %.2g cycles', wc));
    set(gca, 'XTick', 1:numel(symmetryNames), 'XTickLabel', symmetryNames, ...
        'YTick', 1:numel(excitationNames), 'YTickLabel', excitationNames, ...
        'TickLabelInterpreter', 'none');
    nexttile;
    imagesc(log10(condition));
    colorbar;
    title('log10 median condition');
    set(gca, 'XTick', 1:numel(symmetryNames), 'XTickLabel', symmetryNames, ...
        'YTick', 1:numel(excitationNames), 'YTickLabel', excitationNames, ...
        'TickLabelInterpreter', 'none');
    exportgraphics(fig, fullfile(outDir, sprintf('%s_w%.2g.png', prefix, wc)), 'Resolution', 160);
    close(fig);
end
end

function value = percentile(x, p)
x = sort(x(isfinite(x)));
if isempty(x)
    value = NaN;
    return;
end
pos = 1 + (numel(x) - 1) * p / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    value = x(lo);
else
    value = x(lo) + (x(hi) - x(lo)) * (pos - lo);
end
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end
