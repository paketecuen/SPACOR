function data = sp1_synthetic_switched_like_case(config, opts)
%SP1_SYNTHETIC_SWITCHED_LIKE_CASE Generic switched X-like port case.
%
% This generator is intentionally mechanism-agnostic.  It does not assert
% diode/triac/thyristor identity; it only creates OFF intervals where the
% branch is not observable and ON intervals where a canonical physical
% equivalent is visible from the port.

if nargin < 1 || isempty(config)
    config = sp1_default_config();
end
if nargin < 2
    opts = struct();
end

fs = nv3_get_option(opts, 'fs', config.fs);
f0 = nv3_get_option(opts, 'f0', config.f0);
duration = nv3_get_option(opts, 'duration', config.duration);
t = (0:round(duration * fs) - 1).' / fs;

onModel = char(nv3_get_option(opts, 'onModel', 'RL'));
gateKind = char(nv3_get_option(opts, 'gateKind', 'positive_halfcycle'));
truth = complete_truth(nv3_get_option(opts, 'truth', default_truth(onModel)));

[onMask, segmentId] = switched_mask(t, f0, gateKind);
[i, di, qi] = on_current_trace(t, onMask, segmentId, truth);
v = off_voltage(t, f0, truth);
v(onMask) = on_voltage(onModel, truth, i(onMask), di(onMask), qi(onMask));

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = ['switched_' lower(onModel) '_like_' gateKind];
data.label = ['Single-phase switched ' upper(onModel) '-like | ' gateKind];
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'i', i, 'di', di, 'qi', qi, ...
    'onMask', onMask, 'segmentId', segmentId);
data.truth = truth;
data.truth.law = ['switched_' upper(onModel) '_like'];
data.truth.onModel = upper(onModel);
data.truth.gateKind = gateKind;
data.noiseSnrDb = noiseSnrDb;
data.power = struct('measured', vN(:) .* iN(:), 'clean', v(:) .* i(:));

end

function truth = default_truth(onModel)
truth = struct();
truth.V0 = 0.20;
truth.R = 8.0;
truth.L = 25e-3;
truth.C = 150e-6;
truth.Gamma = 1 / truth.C;
truth.currentPeak = 12;
truth.voltageRms = 120;
switch upper(onModel)
    case 'R'
        truth.V0 = 0.05;
    case 'RL'
        truth.V0 = 0.20;
    case 'RC'
        truth.R = 10;
        truth.C = 220e-6;
        truth.Gamma = 1 / truth.C;
    case 'RLC'
        truth.R = 7;
        truth.L = 18e-3;
        truth.C = 180e-6;
        truth.Gamma = 1 / truth.C;
    otherwise
        error('Unknown switched ON model "%s".', onModel);
end
end

function truth = complete_truth(truth)
if ~isfield(truth, 'C') && isfield(truth, 'Gamma')
    truth.C = 1 / truth.Gamma;
end
if ~isfield(truth, 'Gamma') && isfield(truth, 'C')
    truth.Gamma = 1 / truth.C;
end
if ~isfield(truth, 'V0')
    truth.V0 = 0;
end
if ~isfield(truth, 'currentPeak')
    truth.currentPeak = 12;
end
if ~isfield(truth, 'voltageRms')
    truth.voltageRms = 120;
end
end

function [mask, segmentId] = switched_mask(t, f0, gateKind)
phase = mod(f0 * t, 1);
switch gateKind
    case 'positive_halfcycle'
        mask = phase < 0.50;
    case 'late_firing'
        mask = phase >= 0.18 & phase < 0.50;
    case 'burst'
        cycle = floor(f0 * t);
        mask = phase < 0.45 & mod(cycle, 3) ~= 2;
    case 'fast_pwm'
        carrierPerCycle = 12;
        carrierPhase = mod(carrierPerCycle * phase, 1);
        lowFrequencyEnable = phase < 0.90;
        mask = lowFrequencyEnable & carrierPhase < 0.35;
    otherwise
        error('Unknown gate kind "%s".', gateKind);
end
mask = mask(:);
segmentId = zeros(size(mask));
starts = find(diff([false; mask; false]) == 1);
ends = find(diff([false; mask; false]) == -1) - 1;
for k = 1:numel(starts)
    segmentId(starts(k):ends(k)) = k;
end
end

function [i, di, qi] = on_current_trace(t, onMask, segmentId, truth)
i = zeros(size(t));
di = zeros(size(t));
qi = zeros(size(t));
for seg = 1:max(segmentId)
    idx = find(segmentId == seg);
    if numel(idx) < 3
        continue;
    end
    local = linspace(0, 1, numel(idx)).';
    shape = sin(pi * local) + 0.25 * sin(3 * pi * local);
    shape = shape / max(max(abs(shape)), eps);
    i(idx) = truth.currentPeak * shape;
    di(idx) = gradient(i(idx), t(idx));
    qi(idx) = cumtrapz(t(idx), i(idx));
end
onMask = onMask(:);
i(~onMask) = 0;
di(~onMask) = 0;
qi(~onMask) = 0;
end

function v = on_voltage(onModel, truth, i, di, qi)
v = truth.V0 + truth.R * i;
switch upper(onModel)
    case 'R'
        return;
    case 'RL'
        v = v + truth.L * di;
    case 'RC'
        v = v + truth.Gamma * qi;
    case 'RLC'
        v = v + truth.L * di + truth.Gamma * qi;
    otherwise
        error('Unknown switched ON model "%s".', onModel);
end
end

function v = off_voltage(t, f0, truth)
v = truth.voltageRms * sqrt(2) * cos(2 * pi * f0 * t + 0.10);
end

function varargout = add_noise(seed, snrDb, varargin)
rng(seed);
varargout = cell(size(varargin));
for idx = 1:numel(varargin)
    x = varargin{idx};
    if isfinite(snrDb)
        sigma = sqrt(mean(x(:) .^ 2)) * 10 ^ (-snrDb / 20);
        varargout{idx} = x + sigma * randn(size(x));
    else
        varargout{idx} = x;
    end
end
end
