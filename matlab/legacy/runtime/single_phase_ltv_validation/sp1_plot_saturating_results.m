function fig = sp1_plot_saturating_results(data, windowResult, fluxResult, outFile, opts)
%SP1_PLOT_SATURATING_RESULTS Plot LTV saturating-inductor diagnostics.

if nargin < 5
    opts = struct();
end
margin = nv3_get_option(opts, 'truthMargin', 0.10);

W = windowResult.windows;
valid = W.Valid > 0;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1600, 1050]);
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
title('Measured waveforms');

nexttile;
plot(data.clean.i, data.clean.Linc, '.', 'MarkerSize', 3);
grid on;
xlabel('i (A)');
ylabel('L_{inc}(i) (H)');
title('True incremental inductance along trajectory');

nexttile;
h = gobjects(0);
hold on;
legendLabels = {};
if any(~valid)
    h(end + 1) = plot(W.CenterTime(~valid), W.L(~valid), '.', ...
        'Color', [0.7 0.7 0.7]);
    legendLabels{end + 1} = 'invalid';
end
h(end + 1) = plot(W.CenterTime(valid), W.L(valid), '-o', ...
    'LineWidth', 1.1, 'MarkerSize', 3);
legendLabels{end + 1} = 'estimated L_w';
h(end + 1) = plot(W.CenterTime, W.TruthLMean, 'k--', 'LineWidth', 1.0);
legendLabels{end + 1} = 'mean true L_{inc}';
grid on;
xlabel('Time (s)');
ylabel('L_w (H)');
title('Windowed local LTI inductance');
legend(h, legendLabels, 'Location', 'best');
apply_zoom([W.L(valid); W.TruthLMean(:)], median(W.TruthLMean, 'omitnan'), margin);

nexttile;
plot(W.CenterTime(valid), W.R(valid), '-o', 'LineWidth', 1.1, 'MarkerSize', 3);
hold on;
yline(data.truth.R, 'k--', 'Truth R', 'LineWidth', 1.0);
grid on;
xlabel('Time (s)');
ylabel('R_w (ohm)');
title('Windowed local resistance');
apply_zoom(W.R(valid), data.truth.R, margin);

nexttile;
plot(data.clean.i, data.clean.lambda, '.', 'Color', [0.75 0.75 0.75], ...
    'MarkerSize', 3);
hold on;
plot(data.i, fluxResult.lambdaFit, '.', 'MarkerSize', 3);
grid on;
xlabel('i (A)');
ylabel('\lambda (Wb-turn)');
title(sprintf('Flux curve: R estimate %.4g ohm', fluxResult.R));
legend({'truth samples', 'estimated samples'}, 'Location', 'best');

nexttile;
plot(fluxResult.curve.Current, fluxResult.curve.LincTruth, 'k--', ...
    'LineWidth', 1.2);
hold on;
plot(fluxResult.curve.Current, fluxResult.curve.LincEstimate, ...
    'LineWidth', 1.2);
grid on;
xlabel('i (A)');
ylabel('d\lambda/di (H)');
title('Constitutive incremental inductance');
legend({'truth', 'estimated'}, 'Location', 'best');

sgtitle(sprintf('%s | SNR %s | valid windows %.1f%% | flux error %.3g', ...
    data.label, snr_label(data.noiseSnrDb), 100 * mean(valid), ...
    fluxResult.lambdaRelativeError));

if nargin >= 4 && ~isempty(outFile)
    exportgraphics(fig, outFile, 'Resolution', 160);
end

end

function apply_zoom(values, center, margin)
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
