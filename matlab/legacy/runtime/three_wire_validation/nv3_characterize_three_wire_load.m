function report = nv3_characterize_three_wire_load(data, config, opts)
%NV3_CHARACTERIZE_THREE_WIRE_LOAD Physical black-box report for 3-wire loads.
%
% This is the Phase 8 public API.  It wraps the adaptive geometric
% identifier and returns an engineering report rather than only raw
% window-level estimates.
%
% Minimum input:
%   data.vab, data.vbc, data.ia, data.ib
%
% Optional:
%   data.vca, data.ic, data.t, data.fs, data.f0, data.label
%
% Alternative voltage input:
%   data.va, data.vb, data.vc

if nargin < 2 || isempty(config)
    config = nv3_default_config();
end
if nargin < 3
    opts = struct();
end

data = normalize_three_wire_data(data, config);
opts = with_defaults(opts, config);

method = nv3_identify_adaptive_geometric(data, config, opts);
input = input_summary(data);
recommendation = make_recommendation(method, opts, input);

report = struct();
report.input = input;
report.recommendation = recommendation;
report.ranking = summarize_candidate_ranking(method.candidates);
report.windows = method.windows;
report.candidates = method.candidates;
report.parameterSummary = method.parameterSummary;
report.method = method;
report.flow = characterization_flow();
report.options = opts;
report.config = config;

end

function data = normalize_three_wire_data(data, config)
if ~isfield(data, 'vab') || ~isfield(data, 'vbc')
    if all(isfield(data, {'va', 'vb', 'vc'}))
        data.vab = data.va(:) - data.vb(:);
        data.vbc = data.vb(:) - data.vc(:);
    else
        error('Data must contain vab/vbc or va/vb/vc.');
    end
end
data.vab = data.vab(:);
data.vbc = data.vbc(:);
if ~isfield(data, 'vca') || isempty(data.vca)
    data.vca = -data.vab - data.vbc;
else
    data.vca = data.vca(:);
end

if ~isfield(data, 'ia') || ~isfield(data, 'ib')
    error('Data must contain ia and ib.  Field ic is optional.');
end
data.ia = data.ia(:);
data.ib = data.ib(:);
if ~isfield(data, 'ic') || isempty(data.ic)
    data.ic = -data.ia - data.ib;
else
    data.ic = data.ic(:);
end

if ~isfield(data, 'fs') || isempty(data.fs)
    data.fs = config.fs;
end
if ~isfield(data, 'f0') || isempty(data.f0)
    data.f0 = config.f0;
end
if ~isfield(data, 't') || isempty(data.t)
    data.t = (0:numel(data.vab) - 1).' / data.fs;
else
    data.t = data.t(:);
end
if ~isfield(data, 'label') || isempty(data.label)
    data.label = 'three-wire black-box load';
end
if ~isfield(data, 'model') || isempty(data.model)
    data.model = 'three_wire_black_box';
end

sizes = [numel(data.vab), numel(data.vbc), numel(data.vca), ...
    numel(data.ia), numel(data.ib), numel(data.ic), numel(data.t)];
if any(sizes ~= sizes(1))
    error('Voltage, current, and time fields must have the same length.');
end

if ~isfield(data, 'power') || ~isfield(data.power, 'measured')
    data.power = struct();
    data.power.measured = three_wire_power(data.vab, data.vbc, ...
        data.ia, data.ib);
end
end

function opts = with_defaults(opts, config)
opts.candidateMode = nv3_get_option(opts, 'candidateMode', 'mixed_blackbox');
opts.coverageGate = nv3_get_option(opts, 'coverageGate', 0.80);
opts.dominantFractionGate = nv3_get_option(opts, 'dominantFractionGate', 0.50);
opts.kvlResidualGate = nv3_get_option(opts, 'kvlResidualGate', 5e-3);
opts.kclResidualGate = nv3_get_option(opts, 'kclResidualGate', 5e-3);
opts.outsideLibraryResidualGate = nv3_get_option(opts, ...
    'outsideLibraryResidualGate', 5e-2);
opts.windowCycles = nv3_get_option(opts, 'windowCycles', config.windowCycles);
opts.hopCycles = nv3_get_option(opts, 'hopCycles', config.hopCycles);
opts.adaptiveHarmonics = nv3_get_option(opts, 'adaptiveHarmonics', ...
    config.adaptiveHarmonics);
opts.harmonics = nv3_get_option(opts, 'harmonics', config.harmonics);
opts.residualGate = nv3_get_option(opts, 'residualGate', []);
opts.energyGate = nv3_get_option(opts, 'energyGate', []);
opts.selectionRelativeTolerance = nv3_get_option(opts, ...
    'selectionRelativeTolerance', config_value(config, ...
    'selectionRelativeTolerance', 0.10));
end

function recommendation = make_recommendation(method, opts, input)
selected = method.windows.Selected > 0;
[dominantModel, dominantFraction] = dominant_column(method.windows, ...
    'SelectedModel', selected);
[dominantTopology, ~] = dominant_column(method.windows, ...
    'SelectedTopology', selected);
[dominantTerms, ~] = dominant_column(method.windows, ...
    'SelectedTerms', selected);
[dominantFamily, ~] = dominant_column(method.windows, ...
    'SelectedFamily', selected);

coverage = method.summary.SelectedCoverage;
bestUsablePower = best_usable_power(method.candidates);
measurementResidual = max(input.voltageKvlResidual / opts.kvlResidualGate, ...
    input.currentKclResidual / opts.kclResidualGate);
notes = {};
if coverage <= 0
    if measurementResidual > 1 && bestUsablePower <= opts.outsideLibraryResidualGate
        status = 'measurement_inconsistent';
        notes{end + 1} = ...
            'Terminal constraints suggest measurement offset/gain inconsistency.';
    else
        status = 'outside_library';
        notes{end + 1} = ...
            'No passive geometrically adequate candidate was accepted.';
    end
elseif coverage < opts.coverageGate
    status = 'low_coverage';
    notes{end + 1} = 'Accepted windows are below the coverage gate.';
elseif dominantFraction < opts.dominantFractionGate
    status = 'ambiguous_equivalent_family';
    notes{end + 1} = 'Several terminal-equivalent candidates compete.';
else
    status = 'accepted';
end

if contains(dominantModel, 'balanced')
    notes{end + 1} = ...
        'Balanced terminal behavior may not uniquely distinguish delta and wye internals.';
end
if contains(dominantTerms, 'glc') || contains(dominantTerms, 'rlgamma')
    notes{end + 1} = ...
        'Full three-parameter branches require rich geometry in the window.';
end

recommendation = struct();
recommendation.status = status;
recommendation.model = dominantModel;
recommendation.topology = dominantTopology;
recommendation.terms = dominantTerms;
recommendation.family = dominantFamily;
recommendation.dominantFraction = dominantFraction;
recommendation.selectedCoverage = coverage;
recommendation.adequateCoverage = method.summary.AdequateCoverage;
recommendation.medianSelectedN = method.summary.MedianSelectedN;
recommendation.medianCondition = method.summary.MedianCondition;
recommendation.medianSigmaMin = method.summary.MedianSigmaMin;
recommendation.medianGeometricVolume = method.summary.MedianGeometricVolume;
recommendation.medianEquationResidual = method.summary.MedianEquationResidual;
recommendation.medianPowerResidual = method.summary.MedianPowerResidual;
recommendation.medianEnergyResidual = method.summary.MedianEnergyResidual;
recommendation.residualGate = method.residualGate;
recommendation.energyGate = method.energyGate;
recommendation.bestUsablePowerResidual = bestUsablePower;
recommendation.voltageKvlResidual = input.voltageKvlResidual;
recommendation.currentKclResidual = input.currentKclResidual;
recommendation.notes = notes(:);
end

function ranking = summarize_candidate_ranking(candidates)
if isempty(candidates)
    ranking = table();
    return;
end
names = unique(candidates.CandidateModel, 'stable');
rows = cell(numel(names), 12);
for idx = 1:numel(names)
    name = names{idx};
    mask = strcmp(candidates.CandidateModel, name);
    usable = mask & candidates.Usable > 0;
    rows(idx, :) = {name, first_value(candidates.CandidateTopology(mask)), ...
        first_value(candidates.CandidateTerms(mask)), ...
        first_value(candidates.CandidateFamily(mask)), ...
        first_finite(candidates.NParameters(mask)), mean(usable(mask)), ...
        median_finite(candidates.ConditionNumber(usable)), ...
        median_finite(candidates.SigmaMin(usable)), ...
        median_finite(candidates.EquationResidual(usable)), ...
        median_finite(candidates.PowerResidual(usable)), ...
        median_finite(candidates.EnergyResidual(usable)), ...
        median_finite(candidates.PowerResidual(mask))};
end
ranking = cell2table(rows, 'VariableNames', { ...
    'CandidateModel', 'CandidateTopology', 'CandidateTerms', ...
    'CandidateFamily', 'NParameters', 'UsableCoverage', ...
    'MedianConditionUsable', 'MedianSigmaMinUsable', ...
    'MedianEquationResidualUsable', 'MedianPowerResidualUsable', ...
    'MedianEnergyResidualUsable', 'MedianPowerResidualAll'});
ranking = sortrows(ranking, {'MedianPowerResidualUsable', ...
    'MedianEquationResidualUsable', 'NParameters'});
end

function summary = input_summary(data)
voltageScale = max([rms_value(data.vab), rms_value(data.vbc), ...
    rms_value(data.vca), eps]);
currentScale = max([rms_value(data.ia), rms_value(data.ib), ...
    rms_value(data.ic), eps]);
summary = struct();
summary.label = char(data.label);
summary.samples = numel(data.t);
summary.fs = data.fs;
summary.f0 = data.f0;
summary.duration = data.t(end) - data.t(1);
summary.voltageKvlResidual = rms_value(data.vab + data.vbc + data.vca) ...
    / voltageScale;
summary.currentKclResidual = rms_value(data.ia + data.ib + data.ic) ...
    / currentScale;
summary.powerRms = rms_value(data.power.measured);
end

function flow = characterization_flow()
flow = {
    'Normalize terminal three-wire measurements.'
    'Fit each window with adaptive harmonic coordinates.'
    'Evaluate passive delta-parallel and wye-series candidate dictionaries.'
    'Reject candidates unsupported by rank, condition, passivity, or residuals.'
    'Select the simplest adequate terminal physical equivalent.'
    'Report outside-library when no pure candidate is physically adequate.'
    };
end

function [value, fraction] = dominant_column(tableIn, columnName, mask)
value = 'reject';
fraction = 0;
if nargin < 3 || isempty(mask)
    mask = true(height(tableIn), 1);
end
if ~any(mask)
    return;
end
values = string(tableIn.(columnName)(mask));
values = values(values ~= "");
if isempty(values)
    return;
end
[uniqueValues, ~, groupIdx] = unique(values, 'stable');
counts = accumarray(groupIdx, 1);
[count, idx] = max(counts);
value = char(uniqueValues(idx));
fraction = count / numel(values);
end

function value = first_value(values)
if isempty(values)
    value = '';
else
    value = char(string(values(1)));
end
end

function value = first_finite(values)
values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = values(1);
end
end

function value = median_finite(values)
values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = median(values);
end
end

function value = rms_value(x)
x = x(:);
value = sqrt(mean(x .^ 2, 'omitnan'));
end

function value = best_usable_power(candidates)
if isempty(candidates) || height(candidates) == 0
    value = Inf;
    return;
end
values = candidates.PowerResidual(candidates.Usable > 0);
values = values(isfinite(values));
if isempty(values)
    value = Inf;
else
    value = min(values);
end
end

function value = config_value(config, name, defaultValue)
if isfield(config, name) && ~isempty(config.(name))
    value = config.(name);
else
    value = defaultValue;
end
end

function p = three_wire_power(vab, vbc, ia, ib)
p = vab(:) .* ia(:) + vbc(:) .* (ia(:) + ib(:));
end
