function fig = nv3_plot_blackbox_candidate_parameters(candidateTable, caseIndex, ...
    candidateModel, truth, outFile, opts)
%NV3_PLOT_BLACKBOX_CANDIDATE_PARAMETERS Plot candidate theta trajectories.

if nargin < 6
    opts = struct();
end
if ischar(candidateTable) || isstring(candidateTable)
    candidateTable = readtable(candidateTable);
end

mask = candidateTable.CaseIndex == caseIndex & ...
    strcmp(candidateTable.CandidateModel, candidateModel);
rows = candidateTable(mask, :);
if isempty(rows)
    error('No rows found for case %d and candidate "%s".', caseIndex, candidateModel);
end

theta = parse_theta_column(rows.FullTheta);
topology = rows.CandidateTopology{1};
names = parameter_names(topology, size(theta, 2));
labels = strrep(names, 'Gamma', '\Gamma_');
valid = rows.Usable > 0 & all(isfinite(theta), 2);
margin = nv3_get_option(opts, 'truthMargin', 0.10);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1500, 950]);
cols = 3;
tileRows = ceil(numel(names) / cols);
tiledlayout(fig, tileRows, cols, 'Padding', 'compact', 'TileSpacing', 'compact');
for idx = 1:numel(names)
    nexttile;
    values = theta(:, idx);
    plot(rows.CenterTime(valid), values(valid), '-o', ...
        'LineWidth', 1.0, 'MarkerSize', 2.5);
    hold on;
    if isfield(truth, names{idx})
        yline(truth.(names{idx}), '--', 'Truth', 'LineWidth', 0.9);
        apply_truth_zoom(values(valid), truth.(names{idx}), margin);
    end
    grid on;
    xlabel('Time (s)');
    ylabel(labels{idx});
    title(labels{idx});
end
sgtitle(sprintf('Case %d | %s | usable %.1f%%', caseIndex, candidateModel, ...
    100 * mean(valid)));

if nargin >= 5 && ~isempty(outFile)
    exportgraphics(fig, outFile, 'Resolution', 160);
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

function names = parameter_names(topology, n)
switch topology
    case 'delta_parallel'
        names = {'Gab', 'Gbc', 'Gca', 'Gammaab', 'Gammabc', 'Gammaca', ...
            'Cab', 'Cbc', 'Cca'};
    case 'wye_series'
        names = {'Ra', 'Rb', 'Rc', 'La', 'Lb', 'Lc', ...
            'Gammaa', 'Gammab', 'Gammac'};
    otherwise
        names = arrayfun(@(idx) sprintf('theta%d', idx), 1:n, ...
            'UniformOutput', false);
end
names = names(1:min(n, numel(names)));
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
