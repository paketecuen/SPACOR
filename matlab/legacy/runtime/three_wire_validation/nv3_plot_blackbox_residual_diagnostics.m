function fig = nv3_plot_blackbox_residual_diagnostics(candidateTable, ...
    selectionTable, caseIndex, outFile)
%NV3_PLOT_BLACKBOX_RESIDUAL_DIAGNOSTICS Plot best residuals by topology.

if ischar(candidateTable) || isstring(candidateTable)
    candidateTable = readtable(candidateTable);
end
if ischar(selectionTable) || isstring(selectionTable)
    selectionTable = readtable(selectionTable);
end

C = candidateTable(candidateTable.CaseIndex == caseIndex, :);
S = selectionTable(selectionTable.CaseIndex == caseIndex, :);
if isempty(C) || isempty(S)
    error('No rows found for case %d.', caseIndex);
end

topologies = unique(C.CandidateTopology, 'stable');
topologies = topologies(~strcmp(topologies, 'reject'));
t = S.CenterTime;
windows = S.Window;
powerByTopology = NaN(numel(windows), numel(topologies));
equationByTopology = NaN(numel(windows), numel(topologies));
energyByTopology = NaN(numel(windows), numel(topologies));

for topoIdx = 1:numel(topologies)
    topology = topologies{topoIdx};
    for wIdx = 1:numel(windows)
        mask = C.Window == windows(wIdx) & strcmp(C.CandidateTopology, topology) & ...
            C.Usable > 0;
        rows = C(mask, :);
        if isempty(rows)
            continue;
        end
        powerByTopology(wIdx, topoIdx) = min(rows.PowerResidual);
        equationByTopology(wIdx, topoIdx) = min(rows.EquationResidual);
        energyByTopology(wIdx, topoIdx) = min(rows.EnergyResidual);
    end
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1450, 880]);
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot_residual_family(t, powerByTopology, topologies);
hold on;
semilogy(t, S.ResidualGate, 'k--', 'LineWidth', 1.0);
grid on;
xlabel('Time (s)');
ylabel('min power residual');
title('Best terminal power residual by topology');
legend([topologies(:); {'gate'}], 'Location', 'best');

nexttile;
plot_residual_family(t, equationByTopology, topologies);
hold on;
semilogy(t, 2 * S.ResidualGate, 'k--', 'LineWidth', 1.0);
grid on;
xlabel('Time (s)');
ylabel('min equation residual');
title('Best equation residual by topology');
legend([topologies(:); {'2x gate'}], 'Location', 'best');

nexttile;
plot_residual_family(t, energyByTopology, topologies);
hold on;
semilogy(t, S.EnergyGate, 'k--', 'LineWidth', 1.0);
grid on;
xlabel('Time (s)');
ylabel('min energy residual');
title('Best energetic residual by topology');
legend([topologies(:); {'gate'}], 'Location', 'best');

nexttile;
plot_candidate_ranking(C);
title('Best usable candidate medians');

caseName = S.TrueCase{1};
degradation = S.Degradation{1};
sgtitle(sprintf('Case %d | %s | %s', caseIndex, caseName, degradation));

if nargin >= 4 && ~isempty(outFile)
    exportgraphics(fig, outFile, 'Resolution', 160);
end

end

function plot_residual_family(t, values, ~)
colors = lines(size(values, 2));
for idx = 1:size(values, 2)
    semilogy(t, values(:, idx), 'LineWidth', 1.1, 'Color', colors(idx, :));
    hold on;
end
ylim_current = ylim;
ylim([max(min(ylim_current), realmin), max(ylim_current)]);
end

function plot_candidate_ranking(C)
models = unique(C.CandidateModel, 'stable');
medianPower = NaN(numel(models), 1);
for idx = 1:numel(models)
    mask = strcmp(C.CandidateModel, models{idx}) & C.Usable > 0 & ...
        isfinite(C.PowerResidual);
    if any(mask)
        medianPower(idx) = median(C.PowerResidual(mask));
    end
end
[sortedValues, order] = sort(medianPower, 'ascend', 'MissingPlacement', 'last');
order = order(isfinite(sortedValues));
if isempty(order)
    text(0.5, 0.5, 'No usable candidates', 'HorizontalAlignment', 'center');
    axis off;
    return;
end
keep = order(1:min(10, numel(order)));
bar(log10(max(medianPower(keep), realmin)));
grid on;
set(gca, 'XTick', 1:numel(keep), 'XTickLabel', models(keep), ...
    'XTickLabelRotation', 60);
ylabel('log10 median power residual');
end
