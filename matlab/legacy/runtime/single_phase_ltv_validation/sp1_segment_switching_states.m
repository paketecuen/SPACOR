function segmentation = sp1_segment_switching_states(data, config, opts)
%SP1_SEGMENT_SWITCHING_STATES Segment OFF/ON/transition intervals.
%
% Engineering OFF means "not electrically observable at the input port", not
% merely "the current sample is exactly zero".  The default detector therefore
% combines current activity and instantaneous power activity, with robust
% noise-floor guards and minimum-duration cleanup.

if nargin < 2 || isempty(config)
    config = sp1_default_config();
end
if nargin < 3
    opts = struct();
end

i = data.i(:);
v = data.v(:);
t = data.t(:);
n = numel(i);
observableMode = char(nv3_get_option(opts, 'observableMode', ...
    'current_and_power'));
currentThresholdFraction = nv3_get_option(opts, 'currentThresholdFraction', 0.02);
currentThresholdAbs = nv3_get_option(opts, 'currentThresholdAbs', []);
currentNoiseFloorAbs = nv3_get_option(opts, 'currentNoiseFloorAbs', []);
currentNoiseSigmaFactor = nv3_get_option(opts, 'currentNoiseSigmaFactor', 6);
powerThresholdFraction = nv3_get_option(opts, 'powerThresholdFraction', 0.01);
powerThresholdAbs = nv3_get_option(opts, 'powerThresholdAbs', []);
powerNoiseFloorAbs = nv3_get_option(opts, 'powerNoiseFloorAbs', []);
powerNoiseSigmaFactor = nv3_get_option(opts, 'powerNoiseSigmaFactor', 6);
smoothSamples = nv3_get_option(opts, 'smoothSamples', ...
    max(1, round(0.00025 * data.fs)));
guardCycles = nv3_get_option(opts, 'guardCycles', 0.01);
minOnCycles = nv3_get_option(opts, 'minOnCycles', 0.02);
minOffCycles = nv3_get_option(opts, 'minOffCycles', 0.02);

currentEnvelope = moving_rms(i, smoothSamples);
powerEnvelope = moving_average(abs(v .* i), smoothSamples);
currentThresholdAbs = observability_threshold(currentEnvelope, ...
    currentThresholdAbs, currentNoiseFloorAbs, currentThresholdFraction, ...
    currentNoiseSigmaFactor);
powerThresholdAbs = observability_threshold(powerEnvelope, ...
    powerThresholdAbs, powerNoiseFloorAbs, powerThresholdFraction, ...
    powerNoiseSigmaFactor);

currentObservable = currentEnvelope >= currentThresholdAbs;
powerObservable = powerEnvelope >= powerThresholdAbs;
active = combine_observability(currentObservable, powerObservable, ...
    observableMode);

guardSamples = max(1, round(guardCycles * data.fs / data.f0));
minOnSamples = max(4, round(minOnCycles * data.fs / data.f0));
minOffSamples = max(4, round(minOffCycles * data.fs / data.f0));
active = remove_short_true_runs(active, minOnSamples);
active = fill_short_false_gaps(active, minOffSamples);
activeSegments = boolean_segments(active);

state = zeros(n, 1);      % 0 off, 1 on, 2 transition
for s = 1:height(activeSegments)
    a = activeSegments.StartIndex(s);
    b = activeSegments.EndIndex(s);
    transitionStart = max(1, a - guardSamples);
    transitionEnd = min(n, b + guardSamples);
    state(transitionStart:min(n, a + guardSamples - 1)) = 2;
    state(max(1, b - guardSamples + 1):transitionEnd) = 2;
    coreStart = a + guardSamples;
    coreEnd = b - guardSamples;
    if coreEnd - coreStart + 1 >= minOnSamples
        state(coreStart:coreEnd) = 1;
    else
        state(a:b) = 2;
    end
end

onSegments = state_segments(state == 1, t, i);
transitionSegments = state_segments(state == 2, t, i);

segmentation = struct();
segmentation.state = state;
segmentation.onSegments = onSegments;
segmentation.transitionSegments = transitionSegments;
segmentation.thresholdAbs = currentThresholdAbs;
segmentation.currentThresholdAbs = currentThresholdAbs;
segmentation.powerThresholdAbs = powerThresholdAbs;
segmentation.guardSamples = guardSamples;
segmentation.minOnSamples = minOnSamples;
segmentation.minOffSamples = minOffSamples;
segmentation.observability = struct( ...
    'Mode', observableMode, ...
    'CurrentEnvelope', currentEnvelope, ...
    'PowerEnvelope', powerEnvelope, ...
    'CurrentObservable', currentObservable, ...
    'PowerObservable', powerObservable, ...
    'ActiveRaw', combine_observability(currentObservable, powerObservable, ...
        observableMode), ...
    'ActiveClean', active);
segmentation.stats = struct( ...
    'OnCoverage', mean(state == 1), ...
    'OffCoverage', mean(state == 0), ...
    'TransitionCoverage', mean(state == 2), ...
    'ActiveCleanCoverage', mean(active), ...
    'ActiveSegmentCount', height(activeSegments), ...
    'OnSegmentCount', height(onSegments), ...
    'TransitionSegmentCount', height(transitionSegments), ...
    'MedianActiveDuration', median_active_duration(activeSegments, t), ...
    'MedianActiveSamples', median_active_samples(activeSegments), ...
    'MedianOnDuration', median_finite(onSegments.Duration), ...
    'MedianOnSamples', median_finite(onSegments.Samples), ...
    'CurrentThresholdAbs', currentThresholdAbs, ...
    'PowerThresholdAbs', powerThresholdAbs, ...
    'SwitchingEventsPerCycle', height(activeSegments) / ...
        max((t(end) - t(1)) * data.f0, eps));
segmentation.config = config;

end

function value = median_active_duration(segments, t)
if height(segments) == 0
    value = NaN;
    return;
end
durations = zeros(height(segments), 1);
for idx = 1:height(segments)
    durations(idx) = t(segments.EndIndex(idx)) - t(segments.StartIndex(idx));
end
value = median_finite(durations);
end

function value = median_active_samples(segments)
if height(segments) == 0
    value = NaN;
else
    value = median_finite(segments.EndIndex - segments.StartIndex + 1);
end
end

function y = moving_rms(x, width)
y = sqrt(moving_average(x(:) .^ 2, width));
end

function y = moving_average(x, width)
width = max(1, round(width));
kernel = ones(width, 1) / width;
y = conv(x(:), kernel, 'same');
end

function threshold = observability_threshold(envelope, explicitThreshold, ...
    explicitNoiseFloor, fraction, sigmaFactor)
if isempty(explicitThreshold)
    if isempty(explicitNoiseFloor)
        noiseFloor = robust_low_activity_floor(envelope);
    else
        noiseFloor = explicitNoiseFloor;
    end
    threshold = max(fraction * robust_peak(envelope), ...
        sigmaFactor * noiseFloor);
else
    threshold = explicitThreshold;
end
threshold = max(threshold, eps);
end

function floorValue = robust_low_activity_floor(x)
x = abs(x(isfinite(x)));
if isempty(x)
    floorValue = eps;
    return;
end
x = sort(x(:));
n = max(3, round(0.20 * numel(x)));
low = x(1:n);
medLow = median(low);
madLow = median(abs(low - medLow));
floorValue = max(medLow + 1.4826 * madLow, eps);
end

function peak = robust_peak(x)
x = abs(x(isfinite(x)));
if isempty(x)
    peak = eps;
else
    peak = percentile(sort(x(:)), 99);
end
end

function active = combine_observability(currentObservable, powerObservable, mode)
switch mode
    case 'current_only'
        active = currentObservable;
    case 'current_or_power'
        active = currentObservable | powerObservable;
    case 'current_and_power'
        active = currentObservable & powerObservable;
    otherwise
        error('Unknown observableMode "%s".', mode);
end
end

function mask = remove_short_true_runs(mask, minSamples)
segments = boolean_segments(mask);
for idx = 1:height(segments)
    a = segments.StartIndex(idx);
    b = segments.EndIndex(idx);
    if b - a + 1 < minSamples
        mask(a:b) = false;
    end
end
end

function mask = fill_short_false_gaps(mask, minSamples)
segments = boolean_segments(~mask);
for idx = 1:height(segments)
    a = segments.StartIndex(idx);
    b = segments.EndIndex(idx);
    if b - a + 1 < minSamples
        mask(a:b) = true;
    end
end
end

function segments = boolean_segments(mask)
idx = find(diff([false; mask(:); false]) ~= 0);
starts = idx(1:2:end);
ends = idx(2:2:end) - 1;
segments = table(starts, ends, 'VariableNames', ...
    {'StartIndex', 'EndIndex'});
end

function segments = state_segments(mask, t, i)
base = boolean_segments(mask);
rows = cell(height(base), 8);
for idx = 1:height(base)
    a = base.StartIndex(idx);
    b = base.EndIndex(idx);
    rows(idx, :) = {idx, a, b, t(a), t(b), t(b) - t(a), ...
        b - a + 1, max(abs(i(a:b)))};
end
if isempty(rows)
    segments = cell2table(cell(0, 8), 'VariableNames', { ...
        'Segment', 'StartIndex', 'EndIndex', 'StartTime', 'EndTime', ...
        'Duration', 'Samples', 'PeakAbsCurrent'});
else
    segments = cell2table(rows, 'VariableNames', { ...
        'Segment', 'StartIndex', 'EndIndex', 'StartTime', 'EndTime', ...
        'Duration', 'Samples', 'PeakAbsCurrent'});
end
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function value = percentile(x, p)
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
