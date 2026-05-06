function segmentation = sp1_segment_energy_levels(data, config, opts)
%SP1_SEGMENT_ENERGY_LEVELS Segment long captures by slow energy-level changes.
%
% This is not an ON/OFF detector.  It is a macro-segmentation tool for real
% laboratory records whose amplitude, load, or operating point changes over
% seconds.  Each segment should be locally closer to an LTI/LTV windowed
% description than the full record.

if nargin < 2 || isempty(config)
    config = sp1_default_config();
end
if nargin < 3
    opts = struct();
end

v = data.v(:);
i = data.i(:);
t = data.t(:);
if ~isfield(data, 'fs') || isempty(data.fs)
    data.fs = config.fs;
end
if ~isfield(data, 'f0') || isempty(data.f0)
    data.f0 = config.f0;
end
n = numel(t);
blockCycles = nv3_get_option(opts, 'blockCycles', 5);
minSegmentCycles = nv3_get_option(opts, 'minSegmentCycles', 10);
relativeJumpThreshold = nv3_get_option(opts, 'relativeJumpThreshold', 0.20);
blockSamples = max(8, round(blockCycles * data.fs / data.f0));
blockCount = floor(n / blockSamples);
if blockCount < 2
    segmentation = one_segment(data);
    return;
end

blocks = block_table(v, i, t, blockSamples, blockCount);
logI = log(max(blocks.Irms, eps));
logP = log(max(abs(blocks.Pmean), eps));
jumpScore = max(abs(diff(logI)), abs(diff(logP)));
jumpThreshold = log(1 + relativeJumpThreshold);
jumpBlocks = find(jumpScore > jumpThreshold);
segmentBlocks = blocks_to_segments(jumpBlocks, blockCount);
segmentBlocks = merge_short_segments(segmentBlocks, ...
    max(1, ceil(minSegmentCycles / blockCycles)));
segments = segment_table(data, segmentBlocks, blockSamples, n);

segmentation = struct();
segmentation.blocks = blocks;
segmentation.segments = segments;
segmentation.jumpScore = jumpScore;
segmentation.jumpThreshold = jumpThreshold;
segmentation.config = struct('blockCycles', blockCycles, ...
    'minSegmentCycles', minSegmentCycles, ...
    'relativeJumpThreshold', relativeJumpThreshold, ...
    'blockSamples', blockSamples);

end

function segmentation = one_segment(data)
rows = {1, 1, numel(data.t), data.t(1), data.t(end), ...
    data.t(end) - data.t(1), numel(data.t), rms_safe(data.v), ...
    rms_safe(data.i), mean(data.v(:) .* data.i(:), 'omitnan'), ...
    power_factor(data.v, data.i)};
segments = cell2table(rows, 'VariableNames', segment_names());
segmentation = struct('blocks', table(), 'segments', segments, ...
    'jumpScore', [], 'jumpThreshold', NaN, 'config', struct());
end

function blocks = block_table(v, i, t, blockSamples, blockCount)
rows = cell(blockCount, 9);
for b = 1:blockCount
    a = (b - 1) * blockSamples + 1;
    z = b * blockSamples;
    rows(b, :) = {b, a, z, t(a), t(z), t(z) - t(a), ...
        rms_safe(v(a:z)), rms_safe(i(a:z)), ...
        mean(v(a:z) .* i(a:z), 'omitnan')};
end
blocks = cell2table(rows, 'VariableNames', {'Block', 'StartIndex', ...
    'EndIndex', 'StartTime', 'EndTime', 'Duration', 'Vrms', 'Irms', ...
    'Pmean'});
end

function segmentBlocks = blocks_to_segments(jumpBlocks, blockCount)
starts = [1; jumpBlocks(:) + 1];
ends = [jumpBlocks(:); blockCount];
segmentBlocks = [starts, ends];
end

function segmentBlocks = merge_short_segments(segmentBlocks, minBlocks)
idx = 1;
while idx <= size(segmentBlocks, 1)
    lengthBlocks = segmentBlocks(idx, 2) - segmentBlocks(idx, 1) + 1;
    if lengthBlocks >= minBlocks || size(segmentBlocks, 1) == 1
        idx = idx + 1;
        continue;
    end
    if idx == 1
        segmentBlocks(2, 1) = segmentBlocks(1, 1);
        segmentBlocks(1, :) = [];
    else
        segmentBlocks(idx - 1, 2) = segmentBlocks(idx, 2);
        segmentBlocks(idx, :) = [];
    end
end
end

function segments = segment_table(data, segmentBlocks, blockSamples, n)
v = data.v(:);
i = data.i(:);
t = data.t(:);
rows = cell(size(segmentBlocks, 1), numel(segment_names()));
for s = 1:size(segmentBlocks, 1)
    a = (segmentBlocks(s, 1) - 1) * blockSamples + 1;
    z = min(n, segmentBlocks(s, 2) * blockSamples);
    rows(s, :) = {s, a, z, t(a), t(z), t(z) - t(a), z - a + 1, ...
        rms_safe(v(a:z)), rms_safe(i(a:z)), ...
        mean(v(a:z) .* i(a:z), 'omitnan'), power_factor(v(a:z), i(a:z))};
end
segments = cell2table(rows, 'VariableNames', segment_names());
end

function names = segment_names()
names = {'Segment', 'StartIndex', 'EndIndex', 'StartTime', 'EndTime', ...
    'Duration', 'Samples', 'Vrms', 'Irms', 'Pmean', 'PowerFactor'};
end

function pf = power_factor(v, i)
apparent = rms_safe(v) * rms_safe(i);
pf = mean(v(:) .* i(:), 'omitnan') / max(apparent, eps);
end

function value = rms_safe(x)
value = sqrt(mean(x(:) .^ 2, 'omitnan'));
value = max(value, eps);
end
