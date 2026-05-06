function fig = nv3_plot_waveforms_and_parameters(result, data, outFile, opts)
%NV3_PLOT_WAVEFORMS_AND_PARAMETERS Plot waveforms plus windowed parameters.

if nargin < 4
    opts = struct();
end

model = result.model;
windows = result.windows;
valid = windows.Valid > 0;
truth = data.truth;
margin = nv3_get_option(opts, 'truthMargin', 0.10);
[names, labels] = parameter_plot_set(model);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1600, 1200]);
tiledlayout(fig, 5, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile([1 2]);
plot(data.t, data.vab, 'LineWidth', 0.9);
hold on;
plot(data.t, data.vbc, 'LineWidth', 0.9);
plot(data.t, data.vca, 'LineWidth', 0.9);
grid on;
xlabel('Time (s)');
ylabel('Line voltage');
legend({'v_{ab}', 'v_{bc}', 'v_{ca}'}, 'Location', 'best');
title('Measured line voltages');

nexttile;
plot(data.t, data.ia, 'LineWidth', 0.9);
hold on;
plot(data.t, data.ib, 'LineWidth', 0.9);
plot(data.t, data.ic, 'LineWidth', 0.9);
grid on;
xlabel('Time (s)');
ylabel('Line current');
legend({'i_a', 'i_b', 'i_c'}, 'Location', 'best');
title('Measured line currents');

tWindow = windows.CenterTime;
for idx = 1:numel(names)
    nexttile;
    name = names{idx};
    values = windows.(name);
    plot(tWindow(~valid), values(~valid), '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 6);
    hold on;
    plot(tWindow(valid), values(valid), '-o', 'LineWidth', 1.1, 'MarkerSize', 3);
    if isfield(truth, name)
        yline(truth.(name), '--', 'Truth', 'LineWidth', 1.0);
        apply_truth_zoom(values(valid), truth.(name), margin);
    end
    grid on;
    xlabel('Window center (s)');
    ylabel(labels{idx});
    title(labels{idx});
end

coverage = 100 * mean(valid);
sgtitle(sprintf('%s | %s | valid %.1f%%', data.label, model, coverage), ...
    'Interpreter', 'none');

if nargin >= 3 && ~isempty(outFile)
    exportgraphics(fig, outFile, 'Resolution', 160);
end
end

function [names, labels] = parameter_plot_set(model)
switch model
    case {'delta_parallel_gl', 'delta_parallel_glc'}
        names = {'Rab', 'Rbc', 'Rca', 'Lab', 'Lbc', 'Lca'};
        labels = {'R_{ab}', 'R_{bc}', 'R_{ca}', 'L_{ab}', 'L_{bc}', 'L_{ca}'};
        if strcmp(model, 'delta_parallel_glc')
            names = [names, {'Cab', 'Cbc', 'Cca'}];
            labels = [labels, {'C_{ab}', 'C_{bc}', 'C_{ca}'}];
        end
    case 'wye_series_rl'
        names = {'Ra', 'Rb', 'Rc', 'La', 'Lb', 'Lc'};
        labels = {'R_a', 'R_b', 'R_c', 'L_a', 'L_b', 'L_c'};
    otherwise
        error('Unknown model "%s".', model);
end
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
