function regime = sp1_switching_regime_report(segmentation, data, ~, opts)
%SP1_SWITCHING_REGIME_REPORT Classify switching time-scale observability.
%
% The goal is not to identify the semiconductor mechanism.  It is to decide
% whether ON/OFF intervals are long enough to estimate physical ON-state
% parameters, or whether the port should be treated as an averaged/PWM-like
% equivalent over a larger window.

if nargin < 4
    opts = struct();
end

stats = segmentation.stats;
minSamples = nv3_get_option(opts, 'minOnSamplesForIdentification', 20);
minDurationCycles = nv3_get_option(opts, 'minOnDurationCyclesForIdentification', 0.04);
fastEventsPerCycle = nv3_get_option(opts, 'fastEventsPerCycle', 4);
highTransitionCoverage = nv3_get_option(opts, 'highTransitionCoverage', 0.50);
mostlyOnCoverage = nv3_get_option(opts, 'mostlyOnCoverage', 0.90);
mostlyOffCoverage = nv3_get_option(opts, 'mostlyOffCoverage', 0.95);
observableCoverageGate = nv3_get_option(opts, 'observableCoverageGate', 0.05);

medianActiveSamples = stats.MedianActiveSamples;
medianOnSamples = stats.MedianOnSamples;
medianActiveCycles = stats.MedianActiveDuration * data.f0;
medianOnCycles = stats.MedianOnDuration * data.f0;
eventsPerCycle = stats.SwitchingEventsPerCycle;

resolvableByActive = isfinite(medianActiveSamples) && ...
    medianActiveSamples >= minSamples && ...
    isfinite(medianActiveCycles) && medianActiveCycles >= minDurationCycles;
resolvableByCore = isfinite(medianOnSamples) && ...
    medianOnSamples >= minSamples && ...
    isfinite(medianOnCycles) && medianOnCycles >= minDurationCycles;

if (stats.ActiveCleanCoverage < observableCoverageGate && ...
        stats.OnCoverage < observableCoverageGate) || ...
        stats.OffCoverage >= mostlyOffCoverage
    label = 'no_observable_on';
    action = 'Do not estimate internal parameters; report OFF/non-observable.';
elseif stats.OnCoverage >= mostlyOnCoverage && stats.OffCoverage <= 1 - mostlyOnCoverage
    label = 'continuous_or_averaged_on';
    action = 'Prefer non-switching LTI/LTV equivalents; ON/OFF segmentation adds little.';
elseif eventsPerCycle >= fastEventsPerCycle || ...
        stats.TransitionCoverage >= highTransitionCoverage || ...
        (resolvableByActive && ~resolvableByCore)
    label = 'fast_switching_pwm_averaged';
    action = ['Use averaged/windowed physical equivalents; do not force ' ...
        'event-by-event ON parameter identification.'];
elseif resolvableByCore && stats.OnCoverage >= observableCoverageGate && ...
        stats.OffCoverage >= observableCoverageGate
    label = 'resolvable_switching';
    action = 'Use switched_X_like ON/OFF identification.';
else
    label = 'marginal_switching';
    action = ['Switching is detected but ON intervals are weak or short; ' ...
        'increase window, sampling rate, or report low confidence.'];
end

regime = struct();
regime.Label = label;
regime.RecommendedAction = action;
regime.EventsPerCycle = eventsPerCycle;
regime.ActiveCleanCoverage = stats.ActiveCleanCoverage;
regime.OnCoverage = stats.OnCoverage;
regime.OffCoverage = stats.OffCoverage;
regime.TransitionCoverage = stats.TransitionCoverage;
regime.MedianActiveSamples = medianActiveSamples;
regime.MedianOnSamples = medianOnSamples;
regime.MedianActiveCycles = medianActiveCycles;
regime.MedianOnCycles = medianOnCycles;
regime.ResolvableByActiveSegments = resolvableByActive;
regime.ResolvableByOnCore = resolvableByCore;
regime.Thresholds = struct('MinSamples', minSamples, ...
    'MinDurationCycles', minDurationCycles, ...
    'FastEventsPerCycle', fastEventsPerCycle, ...
    'HighTransitionCoverage', highTransitionCoverage);

end
