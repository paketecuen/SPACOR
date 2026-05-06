function data = nv3_apply_measurement_degradation(data, degradation)
%NV3_APPLY_MEASUREMENT_DEGRADATION Apply simple measurement imperfections.
%
% Supported fields:
%   voltageGain      scalar or [vab vbc vca]
%   currentGain      scalar or [ia ib ic]
%   voltageDelaySec  scalar or [vab vbc vca]
%   currentDelaySec  scalar or [ia ib ic]

if nargin < 2 || isempty(degradation)
    return;
end

t = data.t(:);
data.vab = degrade_channel(data.vab, t, gain_at(degradation, 'voltageGain', 1, 1), ...
    gain_at(degradation, 'voltageDelaySec', 0, 1));
data.vbc = degrade_channel(data.vbc, t, gain_at(degradation, 'voltageGain', 1, 2), ...
    gain_at(degradation, 'voltageDelaySec', 0, 2));
data.vca = degrade_channel(data.vca, t, gain_at(degradation, 'voltageGain', 1, 3), ...
    gain_at(degradation, 'voltageDelaySec', 0, 3));
data.ia = degrade_channel(data.ia, t, gain_at(degradation, 'currentGain', 1, 1), ...
    gain_at(degradation, 'currentDelaySec', 0, 1));
data.ib = degrade_channel(data.ib, t, gain_at(degradation, 'currentGain', 1, 2), ...
    gain_at(degradation, 'currentDelaySec', 0, 2));
data.ic = degrade_channel(data.ic, t, gain_at(degradation, 'currentGain', 1, 3), ...
    gain_at(degradation, 'currentDelaySec', 0, 3));

data.power.measured = data.vab(:) .* data.ia(:) + ...
    data.vbc(:) .* (data.ia(:) + data.ib(:));
data.degradation = degradation;

end

function y = degrade_channel(x, t, gain, delaySec)
x = x(:);
if delaySec == 0
    shifted = x;
else
    period = t(end) - t(1) + median(diff(t));
    tq = mod(t - t(1) - delaySec, period) + t(1);
    shifted = interp1(t, x, tq, 'linear', 'extrap');
end
y = gain * shifted;
end

function value = gain_at(s, field, defaultValue, idx)
if ~isstruct(s) || ~isfield(s, field) || isempty(s.(field))
    value = defaultValue;
    return;
end
raw = s.(field);
if numel(raw) == 1
    value = raw;
else
    value = raw(idx);
end
end
