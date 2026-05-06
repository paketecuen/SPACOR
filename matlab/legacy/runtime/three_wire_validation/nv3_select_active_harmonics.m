function selection = nv3_select_active_harmonics(signals, t, f0, candidates, opts)
%NV3_SELECT_ACTIVE_HARMONICS Select harmonics with measurable energy.

if nargin < 5
    opts = struct();
end

signals = double(signals);
t = t(:);
candidates = candidates(:).';
threshold = nv3_get_option(opts, 'activeThreshold', 1e-3);
maxCondition = nv3_get_option(opts, 'maxBasisCondition', 1e10);
selectionMode = char(nv3_get_option(opts, 'harmonicSelectionMode', 'energy'));

amps = zeros(numel(candidates), size(signals, 2));
conditions = zeros(numel(candidates), 1);
omega0 = 2 * pi * f0;

for idx = 1:numel(candidates)
    h = candidates(idx);
    B = [ones(numel(t), 1), cos(h * omega0 * t), sin(h * omega0 * t)];
    s = svd(B, 'econ');
    if isempty(s) || min(s) <= eps(max(s))
        conditions(idx) = Inf;
    else
        conditions(idx) = max(s) / min(s);
    end
    coef = B \ signals;
    amps(idx, :) = hypot(coef(2, :), coef(3, :));
end

if strcmpi(selectionMode, 'residual')
    [selected, selectionTrace] = residual_greedy_selection(signals, t, f0, ...
        candidates, threshold, maxCondition);
else
    energy = max(amps, [], 2);
    scale = max(energy);
    selectionTrace = struct('mode', 'energy', 'residuals', NaN, ...
        'improvements', NaN);
    if ~isfinite(scale) || scale <= eps
        selected = candidates(1);
    else
        active = energy >= threshold * scale & conditions < maxCondition;
        activeCandidates = candidates(active);
        activeEnergy = energy(active);
        [~, order] = sort(activeEnergy, 'descend');
        selected = [];
        for k = 1:numel(order)
            trial = sort([selected, activeCandidates(order(k))]);
            if harmonic_basis_condition(t, f0, trial) < maxCondition
                selected = trial;
            end
        end
        if isempty(selected)
            selected = candidates(1);
        end
    end
end

selection = struct();
selection.selectedHarmonics = selected(:).';
selection.candidateHarmonics = candidates;
selection.amplitudes = amps;
selection.conditions = conditions;
selection.combinedCondition = harmonic_basis_condition(t, f0, selected);
selection.trace = selectionTrace;

end

function [selected, trace] = residual_greedy_selection(signals, t, f0, candidates, threshold, maxCondition)
selected = [];
remaining = candidates(:).';
currentResidual = basis_residual(signals, t, f0, selected);
traceResiduals = currentResidual;
traceImprovements = [];

while ~isempty(remaining)
    bestImprovement = -Inf;
    bestResidual = Inf;
    bestHarmonic = NaN;
    for idx = 1:numel(remaining)
        trial = sort([selected, remaining(idx)]);
        if harmonic_basis_condition(t, f0, trial) >= maxCondition
            continue;
        end
        trialResidual = basis_residual(signals, t, f0, trial);
        improvement = currentResidual - trialResidual;
        if improvement > bestImprovement
            bestImprovement = improvement;
            bestResidual = trialResidual;
            bestHarmonic = remaining(idx);
        end
    end
    if ~isfinite(bestImprovement) || bestImprovement < threshold
        break;
    end
    selected = sort([selected, bestHarmonic]);
    remaining(remaining == bestHarmonic) = [];
    currentResidual = bestResidual;
    traceResiduals(end + 1) = currentResidual; %#ok<AGROW>
    traceImprovements(end + 1) = bestImprovement; %#ok<AGROW>
end

if isempty(selected)
    selected = candidates(1);
end

trace = struct();
trace.mode = 'residual';
trace.residuals = traceResiduals;
trace.improvements = traceImprovements;
end

function residual = basis_residual(signals, t, f0, harmonics)
B = harmonic_basis(t, f0, harmonics);
fit = B * (B \ signals);
residual = norm(signals(:) - fit(:)) / max(norm(signals(:)), eps);
end

function conditionValue = harmonic_basis_condition(t, f0, harmonics)
B = harmonic_basis(t, f0, harmonics);
s = svd(B, 'econ');
if isempty(s) || min(s) <= eps(max(s))
    conditionValue = Inf;
else
    conditionValue = max(s) / min(s);
end
end

function B = harmonic_basis(t, f0, harmonics)
omega0 = 2 * pi * f0;
B = ones(numel(t), 1);
for h = harmonics(:).'
    B = [B, cos(h * omega0 * t), sin(h * omega0 * t)]; %#ok<AGROW>
end
end
