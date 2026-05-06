function fig = sp1_plot_switching_characterization(data, report, outFile, opts)
%SP1_PLOT_SWITCHING_CHARACTERIZATION Plot ON/OFF and parameter estimates.

if nargin < 4
    opts = struct();
end
truthBandFraction = nv3_get_option(opts, 'truthBandFraction', 0.15);

seg = report.segmentation;
est = report.identification.estimates;
valid = est.Valid > 0;
t = data.t(:);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1650, 1100]);
tiledlayout(fig, 4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
yyaxis left;
plot(t, data.i, 'LineWidth', 0.8);
ylabel('i (A)');
yyaxis right;
plot(t, data.v, 'LineWidth', 0.8);
ylabel('v (V)');
grid on;
xlabel('Time (s)');
title('Terminal waveforms');

nexttile;
yyaxis left;
stairs(t, seg.state, 'LineWidth', 1.0);
ylim([-0.2, 2.2]);
set(gca, 'YTick', [0 1 2], 'YTickLabel', {'OFF', 'ON', 'transition'});
ylabel('state');
yyaxis right;
plot(t, seg.observability.CurrentEnvelope ./ ...
    max(seg.currentThresholdAbs, eps), 'LineWidth', 0.8);
hold on;
plot(t, seg.observability.PowerEnvelope ./ ...
    max(seg.powerThresholdAbs, eps), 'LineWidth', 0.8);
yline(1, 'k:', 'LineWidth', 0.9);
ylabel('observable / threshold');
grid on;
xlabel('Time (s)');
title(sprintf('Switching segmentation | ON %.1f%% | OFF %.1f%%', ...
    100 * seg.stats.OnCoverage, 100 * seg.stats.OffCoverage));
legend({'state', 'current activity', 'power activity', 'threshold'}, ...
    'Location', 'best');

nexttile;
hold on;
truthR = truth_value(data, 'R');
draw_truth_band(t, truthR, truthBandFraction, [0.2 0.2 0.2]);
plot_on_segments(est, 'R', ~valid, [0.70 0.70 0.70], ':', 1.0);
plot_on_segments(est, 'R', valid, [0.00 0.45 0.74], '-', 1.6);
draw_truth_line(truthR, 'Truth R');
grid on;
xlim([t(1), t(end)]);
xlabel('Time (s)');
ylabel('R (ohm)');
title('ON-state R estimates');
apply_truth_scale(est.R(valid), truthR, truthBandFraction);

nexttile;
yyaxis left;
hold on;
truthL = truth_value(data, 'L');
draw_truth_band(t, truthL, truthBandFraction, [0.2 0.2 0.2]);
plot_on_segments(est, 'L', ~valid, [0.70 0.70 0.70], ':', 1.0);
plot_on_segments(est, 'L', valid, [0.00 0.45 0.74], '-', 1.6);
draw_truth_line(truthL, 'Truth L');
ylabel('L (H)');
apply_truth_scale(est.L(valid), truthL, truthBandFraction);
yyaxis right;
hold on;
truthV0 = truth_value(data, 'Vd');
draw_truth_band(t, truthV0, truthBandFraction, [0.45 0.18 0.06]);
plot_on_segments(est, 'V0', ~valid, [0.70 0.70 0.70], ':', 1.0);
plot_on_segments(est, 'V0', valid, [0.85 0.33 0.10], '-', 1.4);
draw_truth_line(truthV0, 'Truth Vd');
ylabel('V0 (V)');
apply_truth_scale(est.V0(valid), truthV0, truthBandFraction);
grid on;
xlim([t(1), t(end)]);
xlabel('Time (s)');
title(sprintf('ON-state L and diode drop | residual %.3g', ...
    report.summary.MedianResidual));

sgtitle(sprintf('%s | %s | %s', data.label, report.model, report.reason));

if nargin >= 3 && ~isempty(outFile)
    exportgraphics(fig, outFile, 'Resolution', 160);
end

end

function plot_on_segments(est, field, mask, color, lineStyle, lineWidth)
values = est.(field);
for idx = 1:height(est)
    if mask(idx) && isfinite(values(idx))
        line([est.StartTime(idx), est.EndTime(idx)], ...
            [values(idx), values(idx)], ...
            'Color', color, 'LineStyle', lineStyle, ...
            'LineWidth', lineWidth, 'Marker', 'none');
        plot(mean([est.StartTime(idx), est.EndTime(idx)]), values(idx), ...
            'o', 'Color', color, 'MarkerFaceColor', color, ...
            'MarkerSize', 4, 'HandleVisibility', 'off');
    end
end
end

function draw_truth_band(t, truth, fraction, color)
if ~isfinite(truth)
    return;
end
[lo, hi] = truth_limits(truth, fraction);
patch([t(1), t(end), t(end), t(1)], [lo, lo, hi, hi], color, ...
    'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function draw_truth_line(truth, label)
if isfinite(truth)
    yline(truth, 'k--', label, 'LineWidth', 1.0);
end
end

function apply_truth_scale(values, truth, fraction)
values = values(isfinite(values));
if isempty(values) || ~isfinite(truth)
    return;
end
[lo, hi] = truth_limits(truth, fraction);
ylim([lo, hi]);
end

function [lo, hi] = truth_limits(truth, fraction)
span = fraction * max(abs(truth), eps);
lo = truth - span;
hi = truth + span;
end

function value = truth_value(data, field)
if isfield(data, 'truth') && isfield(data.truth, field)
    value = data.truth.(field);
else
    value = NaN;
end
end
