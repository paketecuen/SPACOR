function result = estimate_derivatives(x, fs, orders, method, opts)
%ESTIMATE_DERIVATIVES Estimate derivatives and return diagnostics.
%
% Supported methods:
%   legacy       - calls historical derivate.m for orders 1..3 when available
%   holoborodko  - Holoborodko first derivative, moment-corrected FIR higher
%   fir_moment   - centered FIR with exact polynomial moments
%   sgolay       - Savitzky-Golay/local polynomial derivative coefficients
%   tikhonov     - finite-difference derivative followed by smoothing
%   harmonic     - harmonic fit and analytic differentiation

if nargin < 4 || isempty(method)
    method = 'fir_moment';
end
if nargin < 5
    opts = struct();
end

x = col(x);
fs = double(fs);
dt = 1 / fs;
orders = unique(round(orders(:))).';
orders = orders(orders >= 0);
maxOrder = max(orders);
method = lower(char(method));

values = cell(maxOrder + 1, 1);
values{1} = x;
diagnostics = empty_diagnostics();

switch method
    case {'provided', 'analytic'}
        provided = get_option(opts, 'provided', {});
        if isempty(provided)
            error('opts.provided is required for method "%s".', method);
        end
        for order = orders
            values{order + 1} = col(provided{order + 1});
            diagnostics = append_diag(diagnostics, order, method, 0, 0, 0, ...
                NaN, NaN, NaN);
        end

    case 'legacy'
        ensure_legacy_derivate_path();
        if exist('derivate', 'file') ~= 2
            error('Historical derivate.m is not on the MATLAB path.');
        end
        for order = orders
            if order == 0
                continue;
            elseif order <= 3
                y = derivate(x, dt, order);
                y = col(y);
                edgeSamples = 12;
            else
                y = x;
                for k = 1:order
                    y = derivate(y, dt, 1);
                end
                y = col(y);
                edgeSamples = 12 * order;
            end
            y = mark_edges(y, edgeSamples);
            values{order + 1} = y;
            diagnostics = append_diag(diagnostics, order, method, NaN, ...
                edgeSamples, NaN, NaN, NaN, NaN);
        end

    case 'tikhonov'
        lambda = get_option(opts, 'lambda', 1e-3);
        allValues = tikhonov_derivatives(x, dt, maxOrder, lambda);
        for order = orders
            values{order + 1} = allValues{order + 1};
            diagnostics = append_diag(diagnostics, order, method, NaN, ...
                2 * order, NaN, lambda, NaN, NaN);
        end

    case 'harmonic'
        f0 = get_option(opts, 'f0', 50);
        [harmonics, harmonicInfo] = resolve_harmonics(x, fs, f0, opts);
        [allValues, fitInfo] = harmonic_derivatives(x, fs, f0, ...
            harmonics, maxOrder);
        harmonicInfo.fit = fitInfo;
        for order = orders
            values{order + 1} = allValues{order + 1};
            diagnostics = append_diag(diagnostics, order, method, NaN, 0, ...
                NaN, NaN, f0, numel(harmonics));
        end

    otherwise
        halfWidth = get_option(opts, 'halfWidth', 15);
        coeffFamily = method;
        coeffOpts = get_option(opts, 'coefficientOptions', struct());
        for order = orders
            if order == 0
                continue;
            end
            [coeffs, info] = derivative_coefficients(order, halfWidth, ...
                coeffFamily, coeffOpts);
            y = conv(x, flipud(coeffs), 'same') / dt^order;
            y = mark_edges(y, info.halfWidth);
            values{order + 1} = y;
            diagnostics = append_diag(diagnostics, order, info.familyUsed, ...
                info.length, info.halfWidth, ...
                info.noiseGainPerDtOrder / dt^order, NaN, NaN, ...
                info.momentError);
        end
end

reference = get_option(opts, 'reference', {});
if ~isempty(reference)
    diagnostics = add_reference_metrics(diagnostics, values, reference, orders);
end

result = struct();
result.method = method;
result.fs = fs;
result.orders = orders;
result.values = values;
result.diagnostics = diagnostics;
if exist('harmonicInfo', 'var')
    result.harmonicInfo = harmonicInfo;
end
result.api = struct( ...
    'name', 'spacor.signal.estimate_derivatives', ...
    'legacyFunction', '', ...
    'legacyFile', '', ...
    'releaseStage', 'public_api');
end

function diagnostics = empty_diagnostics()
diagnostics = table(zeros(0, 1), cell(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'Order', 'Method', 'FilterLength', 'EdgeSamples', ...
    'NoiseGain', 'Lambda', 'F0', 'MomentError', 'NRMSE', 'Gain', ...
    'ValidFraction'});
end

function diagnostics = append_diag(diagnostics, order, method, filterLength, ...
    edgeSamples, noiseGain, lambda, f0, momentError)
row = table(order, {char(method)}, filterLength, edgeSamples, noiseGain, ...
    lambda, f0, momentError, NaN, NaN, NaN, ...
    'VariableNames', diagnostics.Properties.VariableNames);
diagnostics = [diagnostics; row];
end

function y = mark_edges(y, edgeSamples)
y = col(y);
edgeSamples = min(floor(numel(y) / 2), max(0, round(edgeSamples)));
if edgeSamples > 0
    y(1:edgeSamples) = NaN;
    y(end - edgeSamples + 1:end) = NaN;
end
end

function diagnostics = add_reference_metrics(diagnostics, values, reference, orders)
for row = 1:height(diagnostics)
    order = diagnostics.Order(row);
    if ~ismember(order, orders) || order + 1 > numel(reference) || ...
            isempty(reference{order + 1})
        continue;
    end
    y = values{order + 1};
    ref = col(reference{order + 1});
    n = min(numel(y), numel(ref));
    y = y(1:n);
    ref = ref(1:n);
    mask = isfinite(y) & isfinite(ref);
    diagnostics.ValidFraction(row) = mean(mask);
    if any(mask)
        err = y(mask) - ref(mask);
        refRms = sqrt(mean(ref(mask).^2));
        diagnostics.NRMSE(row) = sqrt(mean(err.^2)) / max(refRms, eps);
        diagnostics.Gain(row) = sum(y(mask) .* ref(mask)) / ...
            max(sum(ref(mask).^2), eps);
    end
end
end

function values = tikhonov_derivatives(x, dt, maxOrder, lambda)
n = numel(x);
values = cell(maxOrder + 1, 1);
values{1} = x;

e = ones(n, 1);
D1 = spdiags([-0.5 * e, 0.5 * e], [-1, 1], n, n) / dt;
if n >= 3
    D1(1, 1:3) = [-3, 4, -1] / (2 * dt);
    D1(n, n-2:n) = [1, -4, 3] / (2 * dt);
end
L = spdiags([e, -2 * e, e], -1:1, n, n);
A = speye(n) + lambda * (L' * L);

Dk = speye(n);
for order = 1:maxOrder
    Dk = D1 * Dk;
    fd = Dk * x;
    values{order + 1} = A \ fd;
end
end

function [harmonics, info] = resolve_harmonics(x, fs, f0, opts)
requested = get_option(opts, 'harmonics', [1 3 5 7]);
adaptive = get_option(opts, 'adaptiveHarmonics', false);
if ischar(requested) || isstring(requested)
    if strcmpi(char(requested), 'adaptive')
        adaptive = true;
        requested = get_option(opts, 'harmonicCandidates', ...
            [1 3 5 7 9 11 13]);
    else
        error('Unknown harmonic selector "%s".', char(requested));
    end
end
requested = unique(round(col(requested))).';
requested = requested(requested > 0);
if isempty(requested)
    requested = 1;
end

if adaptive
    candidateHarmonics = get_option(opts, 'harmonicCandidates', requested);
    info = select_active_harmonics(x, fs, f0, candidateHarmonics, opts);
    harmonics = info.selectedHarmonics(:).';
else
    harmonics = requested;
    info = struct('candidateHarmonics', harmonics(:), ...
        'selectedHarmonics', harmonics(:), ...
        'activeMask', true(numel(harmonics), 1), ...
        'adaptive', false);
end
end

function [values, fitInfo] = harmonic_derivatives(x, fs, f0, harmonics, maxOrder)
n = numel(x);
t = (0:n - 1).' / fs;
[cols, basisInfo] = harmonic_basis(n, fs, f0, harmonics);
coef = svd_solve(cols, x);

values = cell(maxOrder + 1, 1);
for order = 0:maxOrder
    y = zeros(n, 1);
    colIdx = 2;
    for h = harmonics(:).'
        omega = 2 * pi * f0 * h;
        theta = omega * t;
        a = coef(colIdx);
        b = coef(colIdx + 1);
        y = y + a * omega^order * cos(theta + order * pi / 2) + ...
            b * omega^order * sin(theta + order * pi / 2);
        colIdx = colIdx + 2;
    end
    if order == 0
        y = y + coef(1);
    end
    values{order + 1} = y;
end

fitInfo = basisInfo;
fitInfo.coefficients = coef;
fitInfo.harmonics = harmonics(:);
fitInfo.reconstructed = cols * coef;
fitInfo.residualRms = sqrt(mean((x - fitInfo.reconstructed) .^ 2));
fitInfo.signalRms = sqrt(mean(x .^ 2));
fitInfo.relativeResidual = fitInfo.residualRms / max(fitInfo.signalRms, eps);
end

function [coeffs, info] = derivative_coefficients(order, halfWidth, family, opts)
if nargin < 3 || isempty(family)
    family = 'fir_moment';
end
if nargin < 4
    opts = struct();
end

order = round(order);
halfWidth = round(halfWidth);
if order < 1
    error('Derivative order must be positive.');
end
if halfWidth < order
    halfWidth = order;
end

family = lower(char(family));
switch family
    case {'fir_moment', 'moment', 'moment_corrected'}
        polyDegree = get_option(opts, 'polyDegree', ...
            min(2 * halfWidth, max(order + 3, 5)));
        coeffs = moment_minimum_norm(order, halfWidth, polyDegree);
        familyUsed = 'fir_moment';
        exactDegree = polyDegree;

    case {'sgolay', 'savitzky_golay', 'savitzky-golay'}
        polyDegree = get_option(opts, 'polyDegree', ...
            min(2 * halfWidth, max(order + 3, 5)));
        coeffs = sgolay_coefficients(order, halfWidth, polyDegree);
        familyUsed = 'sgolay';
        exactDegree = polyDegree;

    case {'holoborodko', 'snrd', 'smooth_noise_robust'}
        if order == 1 && halfWidth >= 2
            coeffs = holoborodko_first(halfWidth);
            familyUsed = 'holoborodko_first';
            exactDegree = 1;
        else
            polyDegree = get_option(opts, 'polyDegree', ...
                min(2 * halfWidth, max(order + 3, 5)));
            coeffs = moment_minimum_norm(order, halfWidth, polyDegree);
            familyUsed = 'fir_moment_higher_order';
            exactDegree = polyDegree;
        end

    otherwise
        error('Unknown derivative coefficient family: %s', family);
end

[coeffs, momentError, momentTable] = normalize_moment(coeffs, order, ...
    halfWidth, exactDegree);

info = struct();
info.familyRequested = family;
info.familyUsed = familyUsed;
info.order = order;
info.halfWidth = halfWidth;
info.length = 2 * halfWidth + 1;
info.exactDegree = exactDegree;
info.offsets = (-halfWidth:halfWidth).';
info.momentError = momentError;
info.momentTable = momentTable;
info.noiseGainPerDtOrder = sqrt(sum(coeffs.^2));
info.l1GainPerDtOrder = sum(abs(coeffs));
end

function coeffs = moment_minimum_norm(order, halfWidth, polyDegree)
offsets = (-halfWidth:halfWidth).';
u = offsets / halfWidth;
polyDegree = min(polyDegree, numel(offsets) - 1);

powers = 0:polyDegree;
A = zeros(numel(powers), numel(offsets));
b = zeros(numel(powers), 1);
for row = 1:numel(powers)
    p = powers(row);
    A(row, :) = u.'.^p;
    if p == order
        b(row) = factorial(order) / halfWidth^order;
    end
end

coeffs = A' * ((A * A') \ b);
end

function coeffs = sgolay_coefficients(order, halfWidth, polyDegree)
offsets = (-halfWidth:halfWidth).';
u = offsets / halfWidth;
polyDegree = min(polyDegree, numel(offsets) - 1);
V = zeros(numel(offsets), polyDegree + 1);
for p = 0:polyDegree
    V(:, p + 1) = u.^p;
end
pinvV = (V' * V) \ V';
coeffs = factorial(order) / halfWidth^order * pinvV(order + 1, :).';
end

function coeffs = holoborodko_first(halfWidth)
M = halfWidth;
m = M - 1;
k = (M:-1:1).';
ck = zeros(size(k));
for idx = 1:numel(k)
    kk = k(idx);
    ck(idx) = (safe_nchoosek(2 * m, m - kk + 1) - ...
        safe_nchoosek(2 * m, m - kk - 1)) / 2^(2 * m + 1);
end
coeffs = -[ck; 0; -flipud(ck)];
end

function value = safe_nchoosek(n, k)
if k < 0 || k > n
    value = 0;
else
    value = exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1));
end
end

function [coeffs, momentError, momentTable] = normalize_moment(coeffs, ...
    order, halfWidth, exactDegree)
coeffs = coeffs(:);
offsets = (-halfWidth:halfWidth).';
target = factorial(order);
moment = sum(coeffs .* offsets.^order);

if moment == 0 || ~isfinite(moment)
    error('Invalid derivative coefficients: order moment is zero or nonfinite.');
end
coeffs = coeffs * (target / moment);

maxPower = min(2 * halfWidth, max([order + 4, 6, exactDegree]));
moments = zeros(maxPower + 1, 1);
targets = zeros(maxPower + 1, 1);
for p = 0:maxPower
    moments(p + 1) = sum(coeffs .* offsets.^p);
    if p == order
        targets(p + 1) = factorial(order);
    end
end
momentError = max(abs(moments(1:exactDegree + 1) - ...
    targets(1:exactDegree + 1)));
momentTable = table((0:maxPower).', moments, targets, moments - targets, ...
    'VariableNames', {'Power', 'Moment', 'Target', 'Error'});
end

function info = select_active_harmonics(signals, fs, f0, candidateHarmonics, opts)
if nargin < 5
    opts = struct();
end

signals = double(signals);
if isvector(signals)
    signals = col(signals);
end
signals = signals(:, any(isfinite(signals), 1));
if isempty(signals)
    signals = zeros(0, 1);
end

fs = double(fs);
f0 = double(f0);
candidateHarmonics = unique(round(col(candidateHarmonics))).';
candidateHarmonics = candidateHarmonics(candidateHarmonics > 0);
candidateHarmonics = candidateHarmonics(f0 * candidateHarmonics < fs / 2);
if isempty(candidateHarmonics)
    candidateHarmonics = 1;
end

activeThreshold = get_option(opts, 'activeThreshold', 1e-3);
referenceMode = lower(char(get_option(opts, 'activeReference', ...
    'fundamental')));
includeFundamental = get_option(opts, 'includeFundamental', true);
maxHarmonics = get_option(opts, 'maxHarmonics', Inf);
maxBasisCondition = get_option(opts, 'maxBasisCondition', ...
    get_option(opts, 'maxHarmonicFitCondition', 1e10));
minSamplesPerColumn = get_option(opts, 'minSamplesPerColumn', 1.05);

n = size(signals, 1);
if n == 0
    info = empty_harmonic_info(candidateHarmonics);
    return;
end

[Phi, basisInfo] = harmonic_basis(n, fs, f0, candidateHarmonics);
coef = svd_solve(Phi, signals);
amplitudesBySignal = harmonic_amplitudes(coef, numel(candidateHarmonics));
amplitude = sqrt(mean(amplitudesBySignal .^ 2, 2));
energy = amplitude .^ 2;

fundIdx = find(candidateHarmonics == 1, 1);
switch referenceMode
    case 'max'
        referenceAmplitude = max(amplitude);
    otherwise
        if isempty(fundIdx)
            referenceAmplitude = max(amplitude);
        else
            referenceAmplitude = amplitude(fundIdx);
            if referenceAmplitude <= max(amplitude) * 1e-6
                referenceAmplitude = max(amplitude);
            end
        end
end
referenceAmplitude = max(referenceAmplitude, eps);
relativeAmplitude = amplitude / referenceAmplitude;
totalEnergy = max(sum(energy), eps);
energyFraction = energy / totalEnergy;

active = relativeAmplitude >= activeThreshold;
if includeFundamental && ~isempty(fundIdx)
    active(fundIdx) = true;
end
if ~any(active)
    [~, strongest] = max(amplitude);
    active(strongest) = true;
end

selected = candidateHarmonics(active);
selectedStrength = amplitude(active);
[selected, selectedStrength] = enforce_harmonic_count(selected, ...
    selectedStrength, maxHarmonics, includeFundamental);
[selected, selectedStrength] = enforce_sample_budget(selected, ...
    selectedStrength, n, minSamplesPerColumn, includeFundamental);
[selected, selectedStrength, selectedCondition, selectedRank] = ...
    enforce_condition(selected, selectedStrength, n, fs, f0, ...
    maxBasisCondition, includeFundamental);

selectedMask = ismember(candidateHarmonics, selected);
selectedPhi = harmonic_basis(n, fs, f0, selected);
selectedCoef = svd_solve(selectedPhi, signals);
reconstructed = selectedPhi * selectedCoef;
residual = signals - reconstructed;
signalRms = sqrt(mean(signals .^ 2, 1)).';
residualRms = sqrt(mean(residual .^ 2, 1)).';

info = struct();
info.candidateHarmonics = candidateHarmonics(:);
info.selectedHarmonics = selected(:);
info.activeMask = selectedMask(:);
info.amplitude = amplitude(:);
info.relativeAmplitude = relativeAmplitude(:);
info.energyFraction = energyFraction(:);
info.referenceAmplitude = referenceAmplitude;
info.fullBasisCondition = basisInfo.conditionNumber;
info.fullBasisRank = basisInfo.rank;
info.selectedBasisCondition = selectedCondition;
info.selectedBasisRank = selectedRank;
info.selectedBasisColumns = 1 + 2 * numel(selected);
info.signalRms = signalRms;
info.residualRms = residualRms;
info.relativeResidual = max(residualRms ./ max(signalRms, eps));
info.options = struct('activeThreshold', activeThreshold, ...
    'activeReference', referenceMode, ...
    'includeFundamental', includeFundamental, ...
    'maxHarmonics', maxHarmonics, ...
    'maxBasisCondition', maxBasisCondition, ...
    'minSamplesPerColumn', minSamplesPerColumn);
end

function info = empty_harmonic_info(candidateHarmonics)
info = struct();
info.candidateHarmonics = candidateHarmonics(:);
info.selectedHarmonics = candidateHarmonics(1);
info.activeMask = false(numel(candidateHarmonics), 1);
info.activeMask(1) = true;
info.amplitude = NaN(numel(candidateHarmonics), 1);
info.relativeAmplitude = NaN(numel(candidateHarmonics), 1);
info.energyFraction = NaN(numel(candidateHarmonics), 1);
info.referenceAmplitude = NaN;
info.fullBasisCondition = Inf;
info.fullBasisRank = 0;
info.selectedBasisCondition = Inf;
info.selectedBasisRank = 0;
info.selectedBasisColumns = NaN;
info.signalRms = NaN;
info.residualRms = NaN;
info.relativeResidual = NaN;
info.options = struct();
end

function [Phi, info] = harmonic_basis(n, fs, f0, harmonics)
t = (0:n - 1).' / fs;
Phi = ones(n, 1);
for h = harmonics(:).'
    theta = 2 * pi * f0 * h * t;
    Phi = [Phi, cos(theta), sin(theta)]; %#ok<AGROW>
end
s = svd(Phi, 'econ');
tol = max(size(Phi)) * eps(max([s; 1]));
rankValue = sum(s > tol);
if isempty(s) || numel(s) < size(Phi, 2) || s(end) <= tol
    conditionNumber = Inf;
else
    conditionNumber = s(1) / s(end);
end
info = struct('singularValues', s, 'rank', rankValue, ...
    'conditionNumber', conditionNumber);
end

function coef = svd_solve(A, B)
[U, S, V] = svd(A, 'econ');
s = diag(S);
tol = max(size(A)) * eps(max([s; 1]));
keep = s > tol;
coef = zeros(size(A, 2), size(B, 2));
if any(keep)
    coef = V(:, keep) * ((U(:, keep)' * B) ./ s(keep));
end
end

function amplitude = harmonic_amplitudes(coef, harmonicCount)
amplitude = zeros(harmonicCount, size(coef, 2));
for idx = 1:harmonicCount
    colIdx = 2 + 2 * (idx - 1);
    amplitude(idx, :) = hypot(coef(colIdx, :), coef(colIdx + 1, :));
end
end

function [selected, strength] = enforce_harmonic_count(selected, strength, ...
    maxHarmonics, includeFundamental)
if ~isfinite(maxHarmonics) || numel(selected) <= maxHarmonics
    return;
end
[~, order] = sort(strength, 'descend');
keep = false(size(selected));
if includeFundamental
    fund = find(selected == 1, 1);
    if ~isempty(fund)
        keep(fund) = true;
    end
end
for idx = order(:).'
    if sum(keep) >= maxHarmonics
        break;
    end
    keep(idx) = true;
end
selected = selected(keep);
strength = strength(keep);
[selected, order] = sort(selected);
strength = strength(order);
end

function [selected, strength] = enforce_sample_budget(selected, strength, ...
    n, minSamplesPerColumn, includeFundamental)
while numel(selected) > 1 && n < ceil(minSamplesPerColumn * ...
        (1 + 2 * numel(selected)))
    [selected, strength] = remove_weakest(selected, strength, ...
        includeFundamental);
end
end

function [selected, strength, conditionNumber, rankValue] = enforce_condition( ...
    selected, strength, n, fs, f0, maxCondition, includeFundamental)
[~, basisInfo] = harmonic_basis(n, fs, f0, selected);
conditionNumber = basisInfo.conditionNumber;
rankValue = basisInfo.rank;
while numel(selected) > 1 && (~isfinite(conditionNumber) || ...
        conditionNumber > maxCondition)
    [selected, strength] = remove_weakest(selected, strength, ...
        includeFundamental);
    [~, basisInfo] = harmonic_basis(n, fs, f0, selected);
    conditionNumber = basisInfo.conditionNumber;
    rankValue = basisInfo.rank;
end
end

function [selected, strength] = remove_weakest(selected, strength, ...
    includeFundamental)
removable = true(size(selected));
if includeFundamental
    removable(selected == 1) = false;
end
if ~any(removable)
    removable = true(size(selected));
end
candidateIdx = find(removable);
[~, local] = min(strength(candidateIdx));
removeIdx = candidateIdx(local);
selected(removeIdx) = [];
strength(removeIdx) = [];
end

function x = col(x)
x = double(x(:));
end

function value = get_option(opts, name, defaultValue)
if nargin < 1 || isempty(opts) || ~isstruct(opts) || ...
        ~isfield(opts, name) || isempty(opts.(name))
    value = defaultValue;
else
    value = opts.(name);
end
end

function ensure_legacy_derivate_path()
if exist('derivate', 'file') == 2
    return;
end
runtimeRoot = fullfile(spacor.core.lab_root(), 'legacy', 'runtime', ...
    'matlab_root_scripts');
if isfile(fullfile(runtimeRoot, 'derivate.m'))
    addpath(runtimeRoot);
    return;
end
matlabRoot = fileparts(spacor.core.lab_root());
if isfile(fullfile(matlabRoot, 'derivate.m'))
    addpath(matlabRoot);
end
end
