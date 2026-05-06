function fig = nv3_plot_characterization_report(report, data, outFile, opts)
%NV3_PLOT_CHARACTERIZATION_REPORT Plot a 3-wire black-box report.
%
% The figure is intentionally engineering-oriented:
%   - terminal waveforms and instantaneous power;
%   - selected physical equivalent by window;
%   - residuals and candidate ranking;
%   - physical parameter trajectories of the dominant accepted model.

if nargin < 4
    opts = struct();
end
if nargin < 3
    outFile = '';
end

data = normalize_plot_data(data, report.config);
truthMargin = nv3_get_option(opts, 'truthMargin', 0.15);
maxRanking = nv3_get_option(opts, 'maxRanking', 10);

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 80, 1800, 1500]);
tiledlayout(fig, 6, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

plot_voltages(data);
plot_currents(data);
plot_power(data);
plot_window_selection(report.windows);
plot_residuals(report.windows);
plot_ranking(report.ranking, maxRanking);
plot_parameters(report, data, truthMargin);

sgtitle(report_title(report), 'Interpreter', 'none');

if ~isempty(outFile)
    outDir = fileparts(outFile);
    if ~isempty(outDir) && ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    exportgraphics(fig, outFile, 'Resolution', 160);
end
end

function data = normalize_plot_data(data, config)
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
if ~isfield(data, 't') || isempty(data.t)
    data.t = (0:numel(data.vab) - 1).' / data.fs;
else
    data.t = data.t(:);
end
if ~isfield(data, 'label') || isempty(data.label)
    data.label = 'three-wire black-box load';
end
end

function plot_voltages(data)
nexttile(1, [1 2]);
plot(data.t, data.vab, 'LineWidth', 0.8);
hold on;
plot(data.t, data.vbc, 'LineWidth', 0.8);
plot(data.t, data.vca, 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Line voltage');
title('Terminal line voltages');
legend({'v_{ab}', 'v_{bc}', 'v_{ca}'}, 'Location', 'best');
end

function plot_currents(data)
nexttile(3);
plot(data.t, data.ia, 'LineWidth', 0.8);
hold on;
plot(data.t, data.ib, 'LineWidth', 0.8);
plot(data.t, data.ic, 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Line current');
title('Terminal line currents');
legend({'i_a', 'i_b', 'i_c'}, 'Location', 'best');
end

function plot_power(data)
nexttile(4, [1 2]);
p = data.vab(:) .* data.ia(:) + data.vbc(:) .* ...
    (data.ia(:) + data.ib(:));
plot(data.t, p, 'k', 'LineWidth', 0.85);
grid on;
xlabel('Time (s)');
ylabel('p(t)');
title('Instantaneous terminal power');
end

function plot_window_selection(windows)
nexttile(6);
models = string(windows.SelectedModel);
models(windows.Selected <= 0) = "reject";
[labels, ~, idx] = unique(models, 'stable');
stairs(windows.CenterTime, idx, 'LineWidth', 1.1);
grid on;
xlabel('Window center (s)');
ylabel('Selected model');
ylim([0.5, max(1.5, numel(labels) + 0.5)]);
set(gca, 'YTick', 1:numel(labels), 'YTickLabel', cellstr(labels));
set(gca, 'TickLabelInterpreter', 'none');
title('Window decision');
end

function plot_residuals(windows)
nexttile(7, [1 2]);
t = windows.CenterTime;
semilogy(t, positive_or_nan(windows.PowerResidual), 'LineWidth', 1.0);
hold on;
semilogy(t, positive_or_nan(windows.EquationResidual), 'LineWidth', 1.0);
semilogy(t, positive_or_nan(windows.EnergyResidual), 'LineWidth', 1.0);
semilogy(t, positive_or_nan(windows.ResidualGate), 'k--', 'LineWidth', 0.9);
semilogy(t, positive_or_nan(windows.EnergyGate), 'k:', 'LineWidth', 0.9);
grid on;
xlabel('Window center (s)');
ylabel('Residual');
title('Accepted-window residual monitors');
legend({'power', 'equation', 'energy', 'residual gate', 'energy gate'}, ...
    'Location', 'best');
end

function plot_ranking(ranking, maxRanking)
nexttile(9);
if isempty(ranking) || height(ranking) == 0
    text(0.5, 0.5, 'No candidates', 'HorizontalAlignment', 'center');
    axis off;
    title('Candidate ranking');
    return;
end
values = ranking.MedianPowerResidualUsable;
valid = isfinite(values);
if ~any(valid)
    text(0.5, 0.5, 'No usable candidates', 'HorizontalAlignment', 'center');
    axis off;
    title('Candidate ranking');
    return;
end
[~, order] = sort(values(valid), 'ascend');
validIdx = find(valid);
keep = validIdx(order(1:min(maxRanking, numel(order))));
bar(log10(max(values(keep), realmin)));
grid on;
set(gca, 'XTick', 1:numel(keep), ...
    'XTickLabel', ranking.CandidateModel(keep), ...
    'XTickLabelRotation', 60);
set(gca, 'TickLabelInterpreter', 'none');
ylabel('log10 median power residual');
title('Best usable candidates');
end

function plot_parameters(report, data, truthMargin)
model = char(report.recommendation.model);
if isempty(model) || strcmp(model, 'reject')
    plot_empty_parameter_tiles('No accepted physical parameters');
    return;
end

windows = report.windows;
candidates = report.candidates;
selectedWindows = windows.Window(windows.Selected > 0 & ...
    strcmp(windows.SelectedModel, model));
rows = candidates(strcmp(candidates.CandidateModel, model), :);
if isempty(rows) || isempty(selectedWindows)
    plot_empty_parameter_tiles('No accepted physical parameters');
    return;
end
mask = rows.Usable > 0 & ismember(rows.Window, selectedWindows);
theta = parse_theta_column(rows.FullTheta);
topology = rows.CandidateTopology{find(mask, 1, 'first')};
termSet = rows.CandidateTerms{find(mask, 1, 'first')};
[names, labels, values, truthScale] = physical_parameter_series(topology, ...
    termSet, theta);
if isempty(names)
    plot_empty_parameter_tiles('No active physical parameters');
    return;
end

time = rows.CenterTime;
for idx = 1:9
    nexttile(9 + idx);
    if idx > numel(names)
        axis off;
        continue;
    end
    y = values(:, idx);
    plot(time(~mask), y(~mask), '.', 'Color', [0.75 0.75 0.75], ...
        'MarkerSize', 5);
    hold on;
    plot(time(mask), y(mask), '-o', 'LineWidth', 1.0, 'MarkerSize', 2.5);
    if isfield(data, 'truth') && isfield(data.truth, names{idx})
        truthValue = truthScale(idx) * data.truth.(names{idx});
        yline(truthValue, '--', 'Truth', 'LineWidth', 0.85);
        apply_truth_zoom(y(mask), truthValue, truthMargin);
    end
    grid on;
    xlabel('Window center (s)');
    ylabel(labels{idx});
    title(labels{idx});
end
end

function plot_empty_parameter_tiles(message)
for idx = 1:9
    nexttile(9 + idx);
    if idx == 5
        text(0.5, 0.5, message, 'HorizontalAlignment', 'center');
    end
    axis off;
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

function [names, labels, values, truthScale] = physical_parameter_series( ...
    topology, termSet, theta)
switch topology
    case 'delta_parallel'
        allNames = {'Rab', 'Rbc', 'Rca', 'Lab', 'Lbc', 'Lca', ...
            'Cab', 'Cbc', 'Cca'};
        allLabels = {'R_{ab}', 'R_{bc}', 'R_{ca}', ...
            'L_{ab}', 'L_{bc}', 'L_{ca}', ...
            'C_{ab} (uF)', 'C_{bc} (uF)', 'C_{ca} (uF)'};
        values = NaN(size(theta, 1), 9);
        values(:, 1:3) = reciprocal(theta(:, 1:3));
        values(:, 4:6) = reciprocal(theta(:, 4:6));
        values(:, 7:9) = 1e6 * theta(:, 7:9);
        allScales = [ones(1, 6), 1e6 * ones(1, 3)];
        active = delta_active_mask(termSet);
    case 'wye_series'
        allNames = {'Ra', 'Rb', 'Rc', 'La', 'Lb', 'Lc', ...
            'Ca', 'Cb', 'Cc'};
        allLabels = {'R_a', 'R_b', 'R_c', ...
            'L_a', 'L_b', 'L_c', 'C_a (uF)', 'C_b (uF)', 'C_c (uF)'};
        values = NaN(size(theta, 1), 9);
        values(:, 1:3) = theta(:, 1:3);
        values(:, 4:6) = theta(:, 4:6);
        values(:, 7:9) = 1e6 * reciprocal(theta(:, 7:9));
        allScales = [ones(1, 6), 1e6 * ones(1, 3)];
        active = wye_active_mask(termSet);
    otherwise
        allNames = {};
        allLabels = {};
        values = [];
        allScales = [];
        active = false(1, 0);
end
names = allNames(active);
labels = allLabels(active);
values = values(:, active);
truthScale = allScales(active);
end

function active = delta_active_mask(termSet)
active = false(1, 9);
if contains(termSet, 'g')
    active(1:3) = true;
end
if contains(termSet, 'l') || contains(termSet, 'gamma')
    active(4:6) = true;
end
if contains(termSet, 'c')
    active(7:9) = true;
end
end

function active = wye_active_mask(termSet)
active = false(1, 9);
if contains(termSet, 'r')
    active(1:3) = true;
end
if contains(termSet, 'l')
    active(4:6) = true;
end
if contains(termSet, 'gamma')
    active(7:9) = true;
end
end

function y = reciprocal(x)
y = NaN(size(x));
mask = isfinite(x) & abs(x) > eps;
y(mask) = 1 ./ x(mask);
end

function values = positive_or_nan(values)
values(~isfinite(values) | values <= 0) = NaN;
end

function apply_truth_zoom(values, truthValue, margin)
values = values(isfinite(values));
if isempty(values) || ~isfinite(truthValue)
    return;
end
p10 = percentile(values, 10);
p90 = percentile(values, 90);
dev = max([abs(p10 - truthValue), abs(p90 - truthValue), ...
    margin * max(abs(truthValue), eps), eps]);
ylim([truthValue - 1.15 * dev, truthValue + 1.15 * dev]);
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

function titleText = report_title(report)
r = report.recommendation;
titleText = sprintf('%s | %s | %s | coverage %.1f%% | p-res %.3g', ...
    report.input.label, r.status, r.model, ...
    100 * r.selectedCoverage, r.medianPowerResidual);
end
