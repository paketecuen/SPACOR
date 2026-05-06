function [dataOut, diagnostics] = sp1_precondition_magnetic_measurements(data, opts)
%SP1_PRECONDITION_MAGNETIC_MEASUREMENTS Explicit offset/drift correction.
%
% Magnetic constitutive models reconstruct flux by integrating v - R i.
% Therefore small measurement offsets create artificial flux ramps.  This
% helper estimates slow measurement components from cycle means, applies the
% correction only when it is plausibly small, and reports the correction.

if nargin < 2
    opts = struct();
end

enabled = nv3_get_option(opts, 'enablePreconditioning', true);
removeOffset = nv3_get_option(opts, 'removeOffset', true);
removeDrift = nv3_get_option(opts, 'removeDrift', false);
maxCorrectionFraction = nv3_get_option(opts, 'maxCorrectionFraction', 0.05);

dataOut = data;
t = data.t(:);
v = data.v(:);
i = data.i(:);
f0 = data.f0;

diagnostics = default_diagnostics();
diagnostics.Enabled = enabled;
diagnostics.RemoveOffset = removeOffset;
diagnostics.RemoveDrift = removeDrift;
diagnostics.MaxCorrectionFraction = maxCorrectionFraction;

if ~enabled || numel(t) < 4 || ~isfinite(f0) || f0 <= 0
    return;
end

degree = double(removeDrift);
if ~removeOffset && ~removeDrift
    return;
end

[tv, vMeans] = cycle_means(t, v, f0);
[ti, iMeans] = cycle_means(t, i, f0);
if numel(tv) < degree + 1 || numel(ti) < degree + 1
    diagnostics.Reason = 'not enough full cycles for preconditioning';
    return;
end

tau = t - mean(t);
vTrend = fitted_trend(tau, tv - mean(t), vMeans, degree);
iTrend = fitted_trend(tau, ti - mean(t), iMeans, degree);

vFraction = max(abs(vTrend)) / rms_safe(v);
iFraction = max(abs(iTrend)) / rms_safe(i);
credible = isfinite(vFraction) && isfinite(iFraction) && ...
    vFraction <= maxCorrectionFraction && iFraction <= maxCorrectionFraction;

diagnostics.VOffset = median(vTrend, 'omitnan');
diagnostics.IOffset = median(iTrend, 'omitnan');
diagnostics.VCorrectionFraction = vFraction;
diagnostics.ICorrectionFraction = iFraction;
diagnostics.Credible = credible;
diagnostics.Reason = sprintf('v correction %.3g rms, i correction %.3g rms', ...
    vFraction, iFraction);

if ~credible
    diagnostics.Applied = false;
    diagnostics.Reason = [diagnostics.Reason ', above credibility gate'];
    return;
end

dataOut.v = v - vTrend;
dataOut.i = i - iTrend;
if isfield(dataOut, 'power') && isstruct(dataOut.power)
    dataOut.power.measured = dataOut.v(:) .* dataOut.i(:);
end
diagnostics.Applied = true;

end

function diagnostics = default_diagnostics()
diagnostics = struct('Enabled', false, 'RemoveOffset', false, ...
    'RemoveDrift', false, 'Applied', false, 'Credible', true, ...
    'MaxCorrectionFraction', NaN, 'VOffset', 0, 'IOffset', 0, ...
    'VCorrectionFraction', 0, 'ICorrectionFraction', 0, ...
    'Reason', 'not requested');
end

function [tc, means] = cycle_means(t, x, f0)
cycle = floor((t - t(1)) * f0);
ids = unique(cycle);
tc = zeros(numel(ids), 1);
means = zeros(numel(ids), 1);
keep = false(numel(ids), 1);
for idx = 1:numel(ids)
    mask = cycle == ids(idx);
    if nnz(mask) < 4
        continue;
    end
    tc(idx) = mean(t(mask), 'omitnan');
    means(idx) = mean(x(mask), 'omitnan');
    keep(idx) = true;
end
tc = tc(keep);
means = means(keep);
end

function trend = fitted_trend(tau, tauCycle, cycleMeans, degree)
if degree == 0
    trend = repmat(median(cycleMeans, 'omitnan'), size(tau));
    return;
end
p = polyfit(tauCycle(:), cycleMeans(:), degree);
trend = polyval(p, tau);
end

function value = rms_safe(x)
value = sqrt(mean(x(:) .^ 2, 'omitnan'));
value = max(value, eps);
end
