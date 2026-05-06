function fig = sp1_plot_bh_parameter_identification(data, result, outFile, opts)
%SP1_PLOT_BH_PARAMETER_IDENTIFICATION Plot B-H parameter estimates.

if nargin < 4
    opts = struct();
end
truthMargin = nv3_get_option(opts, 'truthMargin', 0.10);

t = data.t(:);
lambdaTruth = data.clean.lambda(:) - mean(data.clean.lambda(:), 'omitnan');
lambdaHat = result.lambda(:);
lambdaScale = dot(lambdaHat, lambdaTruth) / max(dot(lambdaHat, lambdaHat), eps);
lambdaHatScaled = lambdaScale * lambdaHat;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1700, 1150]);
tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
yyaxis left;
plot(t, data.i, 'LineWidth', 0.8);
ylabel('i (A)');
yyaxis right;
plot(t, data.v, 'LineWidth', 0.8);
ylabel('v (V)');
grid on;
xlabel('Time (s)');
title('Measured waveforms');

nexttile;
plot(result.ranking.R, result.ranking.Score, 'LineWidth', 1.1);
hold on;
xline(result.R, 'LineWidth', 1.0);
if isfield(data.truth, 'R')
    xline(data.truth.R, 'k--', 'LineWidth', 1.0);
end
grid on;
xlabel('R candidate (ohm)');
ylabel('score');
title('R selection by B-H model fit');
legend({'score', 'selected R', 'truth R'}, 'Location', 'best');

nexttile;
plot(t, lambdaTruth, 'k--', 'LineWidth', 1.0);
hold on;
plot(t, lambdaHatScaled, 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('\lambda (Wb-turn)');
title(sprintf('Reconstructed flux | rel. error %.3g', ...
    result.truth.LambdaRelativeError));
legend({'truth', 'from selected R'}, 'Location', 'best');

nexttile;
plot(data.clean.i, lambdaTruth, '.', 'Color', [0.72 0.72 0.72], ...
    'MarkerSize', 3);
hold on;
plot(data.i, lambdaHatScaled, '.', 'MarkerSize', 3);
grid on;
xlabel('i (A)');
ylabel('\lambda (Wb-turn)');
title('Flux-current loop');
legend({'truth', 'identified'}, 'Location', 'best');

nexttile;
names = {'R', 'L0', 'k3', 'h'};
estimated = [result.R, result.L0, result.saturationCubic, ...
    result.hysteresisCurrent];
truth = [data.truth.R, data.truth.L0, data.truth.saturationCubic, ...
    data.truth.hysteresisCurrent];
bar([truth(:), estimated(:)]);
set(gca, 'XTickLabel', names);
grid on;
ylabel('parameter value');
title('Physical parameter estimates');
legend({'truth', 'estimated'}, 'Location', 'best');
apply_parameter_zoom(truth, estimated, truthMargin);

nexttile;
residual = data.i(:) - result.iFit(:);
plot(t, residual, 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('A');
title(sprintf('B-H current residual | RMS rel. %.3g | loss %.4g J/cycle', ...
    result.residual, result.magneticLossPerCycle));

sgtitle(sprintf('%s | SNR %s | R %.4g | L0 %.4g | k3 %.4g | h %.4g', ...
    data.label, snr_label(data.noiseSnrDb), result.R, result.L0, ...
    result.saturationCubic, result.hysteresisCurrent));

if nargin >= 3 && ~isempty(outFile)
    exportgraphics(fig, outFile, 'Resolution', 160);
end

end

function apply_parameter_zoom(truth, estimated, margin)
values = [truth(:); estimated(:)];
values = values(isfinite(values));
if isempty(values)
    return;
end
lo = min(values);
hi = max(values);
span = hi - lo;
if span <= eps
    span = max(abs(hi), 1);
end
ylim([lo - margin * span, hi + margin * span]);
end

function label = snr_label(snrDb)
if isfinite(snrDb)
    label = sprintf('%.0f dB', snrDb);
else
    label = 'clean';
end
end
