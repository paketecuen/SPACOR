function data = sp1_synthetic_parallel_series_rl_case(config, opts)
%SP1_SYNTHETIC_PARALLEL_SERIES_RL_CASE Synthetic R || (Rs+L) load.
%
% The voltage is prescribed and the current is computed from the exact
% frequency response of the parallel conductance plus the real inductor
% branch:
%
%   Y(jw) = Gp + 1/(Rs + jwL).

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

truth = complete_truth(nv3_get_option(opts, 'truth', default_truth()));
voltageKind = char(nv3_get_option(opts, 'voltageKind', 'multisine'));
voltageComponents = nv3_get_option(opts, 'voltageComponents', ...
    default_voltage_components(voltageKind, truth));

vEval = nv3_harmonic_eval(voltageComponents, t, f0, 1);
v = vEval.values{1};
qv = vEval.primitive;
[i, branchCurrent] = current_from_voltage(voltageComponents, t, f0, truth);

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = ['parallel_series_rl_' voltageKind];
data.label = ['Single-phase R || (Rs+L) | ' voltageKind];
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'dv', vEval.values{2}, 'qv', qv, ...
    'i', i, 'branchCurrent', branchCurrent);
data.truth = truth;
data.truth.law = 'parallel_series_rl';
data.noiseSnrDb = noiseSnrDb;
data.voltageKind = voltageKind;
data.voltageComponents = voltageComponents;
data.harmonics = unique(voltageComponents(:, 2)).';
data.power = struct('measured', vN(:) .* iN(:), 'clean', v(:) .* i(:));

end

function truth = default_truth()
truth = struct();
truth.Rp = 80;
truth.Gp = 1 / truth.Rp;
truth.Rs = 6;
truth.L = 0.080;
truth.voltageRms = 120;
end

function truth = complete_truth(truth)
if ~isfield(truth, 'Gp') && isfield(truth, 'Rp')
    truth.Gp = 1 / truth.Rp;
end
if ~isfield(truth, 'Rp') && isfield(truth, 'Gp')
    truth.Rp = 1 / truth.Gp;
end
end

function components = default_voltage_components(kind, truth)
vpk = truth.voltageRms * sqrt(2);
switch kind
    case 'sine'
        components = [vpk, 1, 0.20];
    case 'multisine'
        components = [ ...
            vpk,        1,  0.20; ...
            0.11 * vpk, 3, -0.35; ...
            0.07 * vpk, 5,  1.05; ...
            0.04 * vpk, 7, -1.45];
    otherwise
        error('Unknown voltage kind "%s".', kind);
end
end

function [i, branchCurrent] = current_from_voltage(components, t, f0, truth)
i = zeros(size(t));
branchCurrent = zeros(size(t));
omega0 = 2 * pi * f0;
for row = 1:size(components, 1)
    amplitude = components(row, 1);
    harmonic = components(row, 2);
    phase = components(row, 3);
    omega = omega0 * harmonic;
    zMag = hypot(truth.Rs, omega * truth.L);
    zPhase = atan2(omega * truth.L, truth.Rs);
    vTerm = amplitude * cos(omega * t + phase);
    branchTerm = amplitude / zMag * cos(omega * t + phase - zPhase);
    branchCurrent = branchCurrent + branchTerm;
    i = i + truth.Gp * vTerm + branchTerm;
end
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
