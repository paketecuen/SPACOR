function fig = nv3_plot_parameter_windows(result, data, outFile, opts)
%NV3_PLOT_PARAMETER_WINDOWS Plot parameter trajectories with truth zoom.

if nargin < 4
    opts = struct();
end

model = result.model;
windows = result.windows;
valid = windows.Valid > 0;
truth = data.truth;
margin = nv3_get_option(opts, 'truthMargin', 0.10);

switch model
    case 'delta_parallel_gl'
        names = {'Rab', 'Rbc', 'Rca', 'Lab', 'Lbc', 'Lca'};
        labels = {'R_{ab}', 'R_{bc}', 'R_{ca}', 'L_{ab}', 'L_{bc}', 'L_{ca}'};
    case 'wye_series_rl'
        names = {'Ra', 'Rb', 'Rc', 'La', 'Lb', 'Lc'};
        labels = {'R_a', 'R_b', 'R_c', 'L_a', 'L_b', 'L_c'};
    otherwise
        error('Unknown model "%s".', model);
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1400, 900]);
t = windows.CenterTime;
for idx = 1:numel(names)
    subplot(3, 2, idx);
    name = names{idx};
    values = windows.(name);
    plot(t(valid), values(valid), '-o', 'LineWidth', 1.2, 'MarkerSize', 3);
    hold on;
    if isfield(truth, name)
        yline(truth.(name), '--', 'Truth', 'LineWidth', 1.0);
        apply_truth_zoom(values(valid), truth.(name), margin);
    end
    grid on;
    xlabel('Time (s)');
    ylabel(labels{idx});
    title(labels{idx});
end
sgtitle(sprintf('%s | %s | valid %.1f%%', data.label, result.model, 100 * mean(valid)));

if nargin >= 3 && ~isempty(outFile)
    exportgraphics(fig, outFile, 'Resolution', 160);
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
