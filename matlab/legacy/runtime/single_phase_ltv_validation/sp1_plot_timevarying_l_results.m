function fig = sp1_plot_timevarying_l_results(data, windowResult, outFile, opts)
%SP1_PLOT_TIMEVARYING_L_RESULTS Plot explicit L(t) identification results.

if nargin < 4
    opts = struct();
end
margin = nv3_get_option(opts, 'truthMargin', 0.10);

W = windowResult.windows;
valid = W.Valid > 0;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1650, 1100]);
tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
yyaxis left;
plot(data.t, data.i, 'LineWidth', 0.9);
ylabel('i (A)');
yyaxis right;
plot(data.t, data.v, 'LineWidth', 0.9);
ylabel('v (V)');
grid on;
xlabel('Time (s)');
title('Terminal waveforms');

nexttile;
yyaxis left;
plot(data.t, data.clean.L, 'LineWidth', 1.0);
ylabel('L(t) (H)');
yyaxis right;
plot(data.t, data.clean.dLdt, 'LineWidth', 1.0);
ylabel('dL/dt (ohm)');
grid on;
xlabel('Time (s)');
title('True LTV law');

nexttile;
h = gobjects(0);
labels = {};
hold on;
if any(~valid)
    h(end + 1) = plot(W.CenterTime(~valid), W.L(~valid), '.', ...
        'Color', [0.72 0.72 0.72]);
    labels{end + 1} = 'invalid';
end
h(end + 1) = plot(W.CenterTime(valid), W.L(valid), '-o', ...
    'LineWidth', 1.1, 'MarkerSize', 3);
labels{end + 1} = 'estimated L_w';
h(end + 1) = plot(W.CenterTime, W.TruthLMean, 'k--', 'LineWidth', 1.0);
labels{end + 1} = 'mean true L(t)';
h(end + 1) = plot(W.CenterTime, W.TruthProjectedL, ':', ...
    'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
labels{end + 1} = 'projected L_{eq}';
grid on;
xlabel('Time (s)');
ylabel('L_w (H)');
title('Windowed inductance');
legend(h, labels, 'Location', 'best');
apply_truth_zoom([W.L(valid); W.TruthLMean; W.TruthProjectedL], ...
    median(W.TruthProjectedL, 'omitnan'), margin);

nexttile;
h = gobjects(0);
labels = {};
hold on;
if any(~valid)
    h(end + 1) = plot(W.CenterTime(~valid), W.RLoss(~valid), '.', ...
        'Color', [0.72 0.72 0.72]);
    labels{end + 1} = 'invalid';
end
h(end + 1) = plot(W.CenterTime(valid), W.RLoss(valid), '-o', ...
    'LineWidth', 1.1, 'MarkerSize', 3);
labels{end + 1} = 'estimated R_w - R';
h(end + 1) = plot(W.CenterTime, W.TruthDldtMean, 'k--', 'LineWidth', 1.0);
labels{end + 1} = 'mean true dL/dt';
h(end + 1) = plot(W.CenterTime, W.TruthProjectedRLoss, ':', ...
    'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
labels{end + 1} = 'projected R_{loss,eq}';
grid on;
xlabel('Time (s)');
ylabel('ohm');
title('Resistance-like LTV term');
legend(h, labels, 'Location', 'best');
apply_truth_zoom([W.RLoss(valid); W.TruthDldtMean; W.TruthProjectedRLoss], ...
    median(W.TruthProjectedRLoss, 'omitnan'), margin);

nexttile;
semilogy(W.CenterTime, W.EquationResidual, '-o', ...
    'LineWidth', 1.0, 'MarkerSize', 3);
hold on;
semilogy(W.CenterTime, W.PowerResidual, '-s', ...
    'LineWidth', 1.0, 'MarkerSize', 3);
grid on;
xlabel('Time (s)');
ylabel('relative residual');
title('Equation and energetic residuals');
legend({'v equation', 'instantaneous power'}, 'Location', 'best');

nexttile;
plot(W.TruthProjectedRLoss(valid), W.RLoss(valid), 'o', ...
    'LineWidth', 1.0, 'MarkerSize', 4);
hold on;
plot_identity(W.TruthProjectedRLoss(valid), W.RLoss(valid));
grid on;
xlabel('projected true R_w - R (ohm)');
ylabel('estimated R_w - R (ohm)');
title('LTV loss-term closure');

sgtitle(sprintf('%s | SNR %s | %.3g cycles | valid %.1f%%', ...
    data.label, snr_label(data.noiseSnrDb), ...
    windowResult.config.windowCycles, 100 * mean(valid)));

if nargin >= 3 && ~isempty(outFile)
    exportgraphics(fig, outFile, 'Resolution', 160);
end

end

function apply_truth_zoom(values, center, margin)
values = values(isfinite(values));
if isempty(values) || ~isfinite(center)
    return;
end
p10 = percentile(values, 10);
p90 = percentile(values, 90);
dev = max([abs(p10 - center), abs(p90 - center), ...
    margin * max(abs(center), eps), eps]);
ylim([center - 1.15 * dev, center + 1.15 * dev]);
end

function plot_identity(x, y)
mask = isfinite(x) & isfinite(y);
x = x(mask);
y = y(mask);
if isempty(x) || isempty(y)
    return;
end
lo = min([x(:); y(:)]);
hi = max([x(:); y(:)]);
if lo == hi
    span = max(abs(lo), 1);
    lo = lo - 0.1 * span;
    hi = hi + 0.1 * span;
end
plot([lo, hi], [lo, hi], 'k--', 'LineWidth', 1.0);
xlim([lo, hi]);
ylim([lo, hi]);
end

function label = snr_label(snrDb)
if isfinite(snrDb)
    label = sprintf('%.0f dB', snrDb);
else
    label = 'clean';
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
