function data = sp1_synthetic_switched_rl_diode_case(config, opts)
%SP1_SYNTHETIC_SWITCHED_RL_DIODE_CASE Ideal diode in series with R+L.
%
% Circuit:
%
%   v(t) -> diode -> R -> L
%
% ON:
%   v = Vd + R i + L di/dt
%
% OFF:
%   i = 0, and R/L are not observable from terminal current.

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
dt = 1 / fs;

truth = nv3_get_option(opts, 'truth', default_truth());
voltageComponents = nv3_get_option(opts, 'voltageComponents', ...
    default_voltage_components(truth));
vEval = nv3_harmonic_eval(voltageComponents, t, f0, 0);
v = vEval.values{1};

[i, on] = simulate_diode_rl(v, dt, truth);

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = 'switched_series_diode_rl';
data.label = 'Single-phase switched diode + series R-L';
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'i', i, 'diodeOn', on);
data.truth = truth;
data.truth.law = 'switched_series_diode_rl';
data.noiseSnrDb = noiseSnrDb;
data.voltageComponents = voltageComponents;
data.harmonics = unique(voltageComponents(:, 2)).';
data.power = struct('measured', vN(:) .* iN(:), 'clean', v(:) .* i(:));

end

function truth = default_truth()
truth = struct();
truth.R = 8.0;
truth.L = 0.025;
truth.Vd = 0.80;
truth.voltageRms = 120;
end

function components = default_voltage_components(truth)
components = [truth.voltageRms * sqrt(2), 1, 0.0];
end

function [i, on] = simulate_diode_rl(v, dt, truth)
n = numel(v);
i = zeros(n, 1);
on = false(n, 1);
conducting = false;
for k = 1:n - 1
    if ~conducting && v(k) > truth.Vd
        conducting = true;
    end

    if conducting
        on(k) = true;
        f1 = diode_rl_rhs(v(k), i(k), truth);
        vMid = 0.5 * (v(k) + v(k + 1));
        f2 = diode_rl_rhs(vMid, i(k) + 0.5 * dt * f1, truth);
        f3 = diode_rl_rhs(vMid, i(k) + 0.5 * dt * f2, truth);
        f4 = diode_rl_rhs(v(k + 1), i(k) + dt * f3, truth);
        iNext = i(k) + dt * (f1 + 2 * f2 + 2 * f3 + f4) / 6;
        if iNext <= 0
            i(k + 1) = 0;
            conducting = false;
        else
            i(k + 1) = iNext;
        end
    else
        i(k + 1) = 0;
    end
end
on(end) = conducting && i(end) > 0;
end

function di = diode_rl_rhs(v, i, truth)
di = (v - truth.Vd - truth.R * i) / truth.L;
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
