function fig = sp1_plot_bh_loop_case(data, outFile, opts)
%SP1_PLOT_BH_LOOP_CASE Plot synthetic B-H loop diagnostics.

if nargin < 3
    opts = struct();
end
R = nv3_get_option(opts, 'R', data.truth.R);

t = data.t(:);
i = data.i(:);
v = data.v(:);
lambdaHat = cumtrapz(t, v - R * i);
lambdaHat = lambdaHat - mean(lambdaHat);
lambdaTruth = data.clean.lambda(:) - mean(data.clean.lambda(:));
scale = dot(lambdaHat, lambdaTruth) / max(dot(lambdaHat, lambdaHat), eps);
lambdaHat = scale * lambdaHat;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1650, 1050]);
tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
yyaxis left;
plot(t, i, 'LineWidth', 0.9);
ylabel('i (A)');
yyaxis right;
plot(t, v, 'LineWidth', 0.9);
ylabel('v (V)');
grid on;
xlabel('Time (s)');
title('Terminal waveforms');

nexttile;
plot(t, lambdaTruth, 'k--', 'LineWidth', 1.1);
hold on;
plot(t, lambdaHat, 'LineWidth', 0.9);
grid on;
xlabel('Time (s)');
ylabel('\lambda (Wb-turn)');
title('Flux from voltage integral');
legend({'truth', 'estimated from v-Ri'}, 'Location', 'best');

nexttile;
plot(data.clean.i, lambdaTruth, 'k--', 'LineWidth', 1.1);
hold on;
plot(i, lambdaHat, '.', 'MarkerSize', 3);
grid on;
xlabel('i (A)');
ylabel('\lambda (Wb-turn)');
title('B-H / flux-current loop');
legend({'truth', 'reconstructed'}, 'Location', 'best');

nexttile;
plot(lambdaTruth, data.clean.iAnhysteretic, 'LineWidth', 1.0);
hold on;
plot(lambdaTruth, data.clean.i, '.', 'MarkerSize', 3);
grid on;
xlabel('\lambda (Wb-turn)');
ylabel('i (A)');
title('Anhysteretic curve and hysteresis branches');
legend({'anhysteretic', 'loop samples'}, 'Location', 'best');

nexttile;
plot(t, data.power.clean, 'LineWidth', 0.8);
hold on;
plot(t, data.power.copper, 'LineWidth', 0.8);
plot(t, data.power.magnetic, 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Power (W)');
title('Power split');
legend({'terminal', 'copper R i^2', 'magnetic i d\lambda/dt'}, ...
    'Location', 'best');

nexttile;
bar([data.energy.MedianMagneticLossPerCycle, ...
    data.energy.MedianCopperLossPerCycle]);
set(gca, 'XTickLabel', {'B-H loop', 'copper'});
grid on;
ylabel('Energy per cycle (J)');
title('Cycle energy');

sgtitle(sprintf('%s | SNR %s | B-H loss %.4g J/cycle', ...
    data.label, snr_label(data.noiseSnrDb), ...
    data.energy.MedianMagneticLossPerCycle));

if nargin >= 2 && ~isempty(outFile)
    exportgraphics(fig, outFile, 'Resolution', 160);
end

end

function label = snr_label(snrDb)
if isfinite(snrDb)
    label = sprintf('%.0f dB', snrDb);
else
    label = 'clean';
end
end
