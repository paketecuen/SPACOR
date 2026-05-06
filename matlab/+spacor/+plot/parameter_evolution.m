function fig = parameter_evolution(report, data, outFile, opts)
%PARAMETER_EVOLUTION Plot windowed parameter trajectories for SPACOR reports.
%
% This is the canonical plotting entry point for quick visual feedback.  It
% supports the current single-phase and three-wire report structures and
% marks invalid windows instead of hiding them.

if nargin < 3
    outFile = '';
end
if nargin < 4
    opts = struct();
end

truthMargin = spacor.core.get_option(opts, 'truthMargin', 0.15);
maxParameters = spacor.core.get_option(opts, 'maxParameters', 9);

if isfield(report, 'windows') && isfield(report, 'ranking')
    plotData = threewire_plot_data(report, data, maxParameters);
    systemLabel = 'three-wire';
else
    plotData = singlephase_plot_data(report, maxParameters);
    systemLabel = 'single-phase';
end

nParams = numel(plotData.parameters);
nRows = max(3, nParams + 2);
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [80, 80, 1550, max(850, 210 * nRows)]);
layout = tiledlayout(fig, nRows, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
title(layout, sprintf('%s | %s | %s', data_label(data), ...
    systemLabel, plotData.selectedModel), 'Interpreter', 'none');

nexttile(layout);
plot_waveforms(data);

for idx = 1:nParams
    nexttile(layout);
    plot_parameter_series(plotData.time, plotData.valid, ...
        plotData.parameters(idx), data, truthMargin);
end

nexttile(layout);
plot_health(plotData);

if ~isempty(outFile)
    outDir = fileparts(outFile);
    if ~isempty(outDir)
        spacor.io.ensure_dir(outDir);
    end
    exportgraphics(fig, outFile, 'Resolution', 170);
end
end

function plotData = singlephase_plot_data(report, maxParameters)
selectedModel = selected_model(report);
[method, selectedModel] = singlephase_method(report, selectedModel);
if isempty(method) || ~isfield(method, 'windows') || isempty(method.windows)
    error('spacor:missingSinglephaseWindows', ...
        'No single-phase window table found for selected model "%s".', ...
        selectedModel);
end
windows = method.windows;
paramNames = singlephase_parameter_names(selectedModel, windows);
paramNames = paramNames(1:min(numel(paramNames), maxParameters));

plotData = struct();
plotData.selectedModel = selectedModel;
plotData.time = windows.CenterTime;
plotData.valid = validity_mask(windows);
plotData.parameters = make_parameters(windows, paramNames);
plotData.condition = table_column_or_nan(windows, 'ConditionNumber');
plotData.equationResidual = table_column_or_nan(windows, 'EquationResidual');
plotData.powerResidual = table_column_or_nan(windows, 'PowerResidual');
end

function plotData = threewire_plot_data(report, data, maxParameters)
selectedModel = selected_model(report);
windows = report.windows;
candidates = report.candidates;

rows = candidates(strcmp(candidates.CandidateModel, selectedModel), :);
if isempty(rows)
    error('spacor:missingThreewireCandidate', ...
        'No candidate rows found for selected model "%s".', selectedModel);
end
usable = rows.Usable > 0;
if any(usable)
    topology = rows.CandidateTopology{find(usable, 1, 'first')};
    termSet = rows.CandidateTerms{find(usable, 1, 'first')};
else
    topology = rows.CandidateTopology{1};
    termSet = rows.CandidateTerms{1};
end

theta = parse_theta_column(rows.FullTheta);
[paramNames, values] = threewire_parameter_series(topology, termSet, theta);
paramNames = paramNames(1:min(numel(paramNames), maxParameters));
values = values(:, 1:numel(paramNames));

parameters = repmat(struct('name', '', 'values', [], 'truthName', '', ...
    'scale', 1, 'label', ''), 0, 1);
for idx = 1:numel(paramNames)
    parameters(idx, 1) = struct('name', paramNames{idx}, ...
        'values', values(:, idx), 'truthName', paramNames{idx}, ...
        'scale', 1, 'label', parameter_label(paramNames{idx}));
end

plotData = struct();
plotData.selectedModel = selectedModel;
plotData.time = rows.CenterTime;
plotData.valid = usable;
plotData.parameters = parameters;
plotData.condition = rows.ConditionNumber;
plotData.equationResidual = rows.EquationResidual;
plotData.powerResidual = rows.PowerResidual;

if isfield(data, 'truth')
    plotData.truth = data.truth;
end
end

function selectedModel = selected_model(report)
if isfield(report, 'canonical') && isfield(report.canonical, 'selectedModel')
    selectedModel = char(report.canonical.selectedModel);
elseif isfield(report, 'recommendation') && isfield(report.recommendation, 'model')
    selectedModel = char(report.recommendation.model);
else
    selectedModel = 'unknown';
end
end

function [method, selectedModel] = singlephase_method(report, selectedModel)
method = [];
if strcmp(selectedModel, 'lti_equivalent_family')
    candidates = report.candidates;
    accepted = strcmp(candidates.Family, 'lti_windowed') & ...
        candidates.Accepted > 0 & isfinite(candidates.ResidualOrScore);
    if any(accepted)
        idx = find(accepted);
        [~, order] = min(candidates.ResidualOrScore(idx));
        selectedModel = candidates.Candidate{idx(order)};
    end
end

fieldName = singlephase_method_field(selectedModel);
if isfield(report, 'methods') && isfield(report.methods, fieldName)
    method = report.methods.(fieldName);
end
end

function fieldName = singlephase_method_field(model)
switch char(model)
    case 'lti_series_rl'
        fieldName = 'seriesRL';
    case 'lti_parallel_rl'
        fieldName = 'parallelRL';
    case 'lti_parallel_series_rl'
        fieldName = 'parallelSeriesRL';
    case 'lti_series_rc'
        fieldName = 'seriesRC';
    case 'lti_parallel_rc'
        fieldName = 'parallelRC';
    case 'lti_capacitor'
        fieldName = 'capacitor';
    case 'lti_series_rlc'
        fieldName = 'seriesRLC';
    case 'lti_parallel_rlc'
        fieldName = 'parallelRLC';
    otherwise
        fieldName = '';
end
end

function names = singlephase_parameter_names(model, windows)
switch char(model)
    case 'lti_series_rl'
        names = {'R', 'L'};
    case 'lti_parallel_rl'
        names = {'G', 'R', 'Gamma', 'L'};
    case 'lti_parallel_series_rl'
        names = {'Gp', 'Rp', 'Rs', 'L'};
    case 'lti_series_rc'
        names = {'R', 'C'};
    case 'lti_parallel_rc'
        names = {'G', 'R', 'C'};
    case 'lti_capacitor'
        names = {'C'};
    case 'lti_series_rlc'
        names = {'R', 'L', 'C'};
    case 'lti_parallel_rlc'
        names = {'G', 'R', 'Gamma', 'L', 'C'};
    otherwise
        names = numeric_parameter_columns(windows);
end
names = names(ismember(names, windows.Properties.VariableNames));
end

function names = numeric_parameter_columns(windows)
blocked = {'Window', 'CenterTime', 'StartTime', 'EndTime', 'Rank', ...
    'ConditionNumber', 'EquationResidual', 'PowerResidual', 'Passive', ...
    'Valid'};
names = {};
for idx = 1:numel(windows.Properties.VariableNames)
    name = windows.Properties.VariableNames{idx};
    if ismember(name, blocked) || startsWith(name, 'Truth') || ...
            startsWith(name, 'Relative')
        continue;
    end
    if isnumeric(windows.(name))
        names{end + 1} = name; %#ok<AGROW>
    end
end
end

function parameters = make_parameters(windows, names)
parameters = repmat(struct('name', '', 'values', [], 'truthName', '', ...
    'scale', 1, 'label', ''), 0, 1);
for idx = 1:numel(names)
    name = names{idx};
    parameters(idx, 1) = struct('name', name, ...
        'values', windows.(name), 'truthName', name, 'scale', 1, ...
        'label', parameter_label(name));
end
end

function valid = validity_mask(windows)
if ismember('Valid', windows.Properties.VariableNames)
    valid = windows.Valid > 0;
elseif ismember('Selected', windows.Properties.VariableNames)
    valid = windows.Selected > 0;
else
    valid = true(height(windows), 1);
end
end

function plot_waveforms(data)
t = time_vector(data);
if isfield(data, 'v') && isfield(data, 'i')
    v = data.v(:);
    i = data.i(:);
    p = v .* i;
    idx = decimated_indices(numel(t), 6000);
    plot(t(idx), v(idx) ./ rms_safe(v), 'LineWidth', 0.85);
    hold on;
    plot(t(idx), i(idx) ./ rms_safe(i), 'LineWidth', 0.85);
    plot(t(idx), p(idx) ./ rms_safe(p), 'LineWidth', 0.85);
    legend({'v/rms(v)', 'i/rms(i)', 'p/rms(p)'}, 'Location', 'best');
    ylabel('normalized');
else
    [vab, vbc, vca, ia, ib, ic] = threewire_signals(data);
    p = vab .* ia + vbc .* (ia + ib);
    idx = decimated_indices(numel(t), 6000);
    plot(t(idx), vab(idx) ./ rms_safe(vab), 'LineWidth', 0.8);
    hold on;
    plot(t(idx), vbc(idx) ./ rms_safe(vbc), 'LineWidth', 0.8);
    plot(t(idx), vca(idx) ./ rms_safe(vca), 'LineWidth', 0.8);
    plot(t(idx), ia(idx) ./ rms_safe(ia), '--', 'LineWidth', 0.8);
    plot(t(idx), ib(idx) ./ rms_safe(ib), '--', 'LineWidth', 0.8);
    plot(t(idx), ic(idx) ./ rms_safe(ic), '--', 'LineWidth', 0.8);
    plot(t(idx), p(idx) ./ rms_safe(p), 'k', 'LineWidth', 0.9);
    legend({'vab', 'vbc', 'vca', 'ia', 'ib', 'ic', 'p'}, ...
        'Location', 'best');
    ylabel('normalized');
end
hold off;
grid on;
xlabel('Time (s)');
title('Terminal waveforms');
end

function plot_parameter_series(time, valid, parameter, data, truthMargin)
y = parameter.values(:);
time = time(:);
valid = valid(:) & isfinite(y);
invalid = ~valid & isfinite(y);
hold on;
truthValue = truth_value(data, parameter.truthName, parameter.scale);
if isfinite(truthValue)
    draw_truth_band(time, truthValue, truthMargin);
end
if any(invalid)
    plot(time(invalid), y(invalid), '.', 'Color', [0.70 0.70 0.70], ...
        'MarkerSize', 8);
end
if any(valid)
    plot(time(valid), y(valid), '-o', 'Color', [0.0 0.45 0.74], ...
        'MarkerFaceColor', [0.0 0.45 0.74], 'LineWidth', 1.1, ...
        'MarkerSize', 3.5);
end

if isfinite(truthValue)
    yline(truthValue, '--', 'Truth', 'LineWidth', 0.9, ...
        'LabelHorizontalAlignment', 'left');
end
medianValue = median_finite(y(valid));
if isfinite(medianValue)
    yline(medianValue, '-', 'Median', 'Color', [0.15 0.15 0.15], ...
        'LineWidth', 0.8, 'LabelHorizontalAlignment', 'right');
end
hold off;
grid on;
xlabel('Window center (s)');
ylabel(parameter.label);
title(parameter.label, 'Interpreter', 'none');
ylim(zoom_limits(y(valid), truthValue, truthMargin));
end

function plot_health(plotData)
t = plotData.time(:);
hasResidual = any(isfinite(plotData.powerResidual)) || ...
    any(isfinite(plotData.equationResidual));
if hasResidual
    semilogy(t, positive_or_nan(plotData.powerResidual), '-o', ...
        'LineWidth', 1.0, 'MarkerSize', 3);
    hold on;
    semilogy(t, positive_or_nan(plotData.equationResidual), '-s', ...
        'LineWidth', 1.0, 'MarkerSize', 3);
else
    plot(t, double(plotData.valid), '-o', 'LineWidth', 1.0);
    hold on;
end
if any(isfinite(plotData.condition))
    yyaxis right;
    semilogy(t, positive_or_nan(plotData.condition), '-.', ...
        'LineWidth', 1.0);
    ylabel('condition');
    yyaxis left;
end
hold off;
grid on;
xlabel('Window center (s)');
ylabel('residual');
title('Window health');
legend({'power residual', 'equation residual', 'condition'}, ...
    'Location', 'best');
end

function [names, values] = threewire_parameter_series(topology, termSet, theta)
switch char(topology)
    case 'delta_parallel'
        allNames = {'Rab', 'Rbc', 'Rca', 'Lab', 'Lbc', 'Lca', ...
            'Cab', 'Cbc', 'Cca'};
        values = NaN(size(theta, 1), 9);
        values(:, 1:3) = reciprocal(theta(:, 1:3));
        values(:, 4:6) = reciprocal(theta(:, 4:6));
        values(:, 7:9) = theta(:, 7:9);
        active = delta_active_mask(termSet);
    case 'wye_series'
        allNames = {'Ra', 'Rb', 'Rc', 'La', 'Lb', 'Lc', ...
            'Ca', 'Cb', 'Cc'};
        values = NaN(size(theta, 1), 9);
        values(:, 1:3) = theta(:, 1:3);
        values(:, 4:6) = theta(:, 4:6);
        values(:, 7:9) = reciprocal(theta(:, 7:9));
        active = wye_active_mask(termSet);
    otherwise
        names = {};
        values = [];
        return;
end
names = allNames(active);
values = values(:, active);
end

function active = delta_active_mask(termSet)
switch char(termSet)
    case 'g'
        active = [true true true false false false false false false];
    case 'gl'
        active = [true true true true true true false false false];
    otherwise
        active = true(1, 9);
end
end

function active = wye_active_mask(termSet)
switch char(termSet)
    case 'r'
        active = [true true true false false false false false false];
    case 'rl'
        active = [true true true true true true false false false];
    otherwise
        active = true(1, 9);
end
end

function theta = parse_theta_column(values)
count = numel(values);
thetaCells = cell(count, 1);
maxLen = 0;
for idx = 1:count
    text = values{idx};
    if iscell(text)
        text = text{1};
    end
    text = strrep(strrep(char(text), '[', ''), ']', '');
    thetaCells{idx} = sscanf(text, '%f').';
    maxLen = max(maxLen, numel(thetaCells{idx}));
end
theta = NaN(count, maxLen);
for idx = 1:count
    theta(idx, 1:numel(thetaCells{idx})) = thetaCells{idx};
end
end

function y = reciprocal(x)
y = NaN(size(x));
mask = isfinite(x) & abs(x) > eps;
y(mask) = 1 ./ x(mask);
end

function value = truth_value(data, name, scale)
value = NaN;
if isfield(data, 'truth') && isfield(data.truth, name)
    value = scale * data.truth.(name);
elseif isfield(data, 'truth')
    alias = truth_alias(name, data.truth);
    if ~isempty(alias)
        value = scale * data.truth.(alias);
    end
end
end

function alias = truth_alias(name, truth)
alias = '';
switch char(name)
    case 'Rs'
        candidates = {'Rser', 'R_s'};
    case 'L'
        candidates = {'Lser', 'L_s'};
    case 'Rp'
        candidates = {'Rpar', 'R_p'};
    otherwise
        candidates = {};
end
for idx = 1:numel(candidates)
    if isfield(truth, candidates{idx})
        alias = candidates{idx};
        return;
    end
end
end

function draw_truth_band(time, truthValue, margin)
if isempty(time) || ~isfinite(truthValue)
    return;
end
lo = truthValue * (1 - margin);
hi = truthValue * (1 + margin);
if truthValue < 0
    lo = truthValue * (1 + margin);
    hi = truthValue * (1 - margin);
end
x = [min(time), max(time), max(time), min(time)];
y = [lo, lo, hi, hi];
patch(x, y, [0.85 0.90 1.0], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.45);
end

function limits = zoom_limits(values, truthValue, margin)
values = values(isfinite(values));
if isfinite(truthValue)
    lo = truthValue * (1 - margin);
    hi = truthValue * (1 + margin);
    if truthValue < 0
        lo = truthValue * (1 + margin);
        hi = truthValue * (1 - margin);
    end
    if ~isempty(values)
        p10 = percentile(values, 10);
        p90 = percentile(values, 90);
        lo = min(lo, p10);
        hi = max(hi, p90);
    end
else
    if isempty(values)
        limits = [0 1];
        return;
    end
    lo = percentile(values, 5);
    hi = percentile(values, 95);
end
if ~isfinite(lo) || ~isfinite(hi) || lo == hi
    center = median_finite([values(:); truthValue]);
    span = max(abs(center), 1);
    limits = [center - 0.1 * span, center + 0.1 * span];
    return;
end
pad = 0.12 * max(abs(hi - lo), eps);
limits = [lo - pad, hi + pad];
end

function label = parameter_label(name)
switch char(name)
    case {'R', 'Rp', 'Rs', 'Ra', 'Rb', 'Rc', 'Rab', 'Rbc', 'Rca'}
        label = [char(name) ' (ohm)'];
    case {'L', 'La', 'Lb', 'Lc', 'Lab', 'Lbc', 'Lca', 'Lser'}
        label = [char(name) ' (H)'];
    case {'C', 'Ca', 'Cb', 'Cc', 'Cab', 'Cbc', 'Cca'}
        label = [char(name) ' (F)'];
    case {'G', 'Gp'}
        label = [char(name) ' (S)'];
    case 'Gamma'
        label = 'Gamma (1/H)';
    otherwise
        label = char(name);
end
end

function [vab, vbc, vca, ia, ib, ic] = threewire_signals(data)
vab = data.vab(:);
vbc = data.vbc(:);
if isfield(data, 'vca') && ~isempty(data.vca)
    vca = data.vca(:);
else
    vca = -vab - vbc;
end
ia = data.ia(:);
ib = data.ib(:);
if isfield(data, 'ic') && ~isempty(data.ic)
    ic = data.ic(:);
else
    ic = -ia - ib;
end
end

function t = time_vector(data)
if isfield(data, 't') && ~isempty(data.t)
    t = data.t(:);
else
    fs = data.fs;
    if isfield(data, 'v')
        n = numel(data.v);
    else
        n = numel(data.vab);
    end
    t = (0:n - 1).' / fs;
end
end

function idx = decimated_indices(n, maxPoints)
if n <= maxPoints
    idx = (1:n).';
else
    idx = unique(round(linspace(1, n, maxPoints))).';
end
end

function value = rms_safe(x)
value = sqrt(mean(x(:) .^ 2, 'omitnan'));
value = max(value, eps);
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function y = table_column_or_nan(tableIn, name)
if ismember(name, tableIn.Properties.VariableNames)
    y = tableIn.(name);
else
    y = NaN(height(tableIn), 1);
end
end

function y = positive_or_nan(x)
y = x;
y(~isfinite(y) | y <= 0) = NaN;
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

function label = data_label(data)
if isfield(data, 'label') && ~isempty(data.label)
    label = char(data.label);
else
    label = 'SPACOR load';
end
end
