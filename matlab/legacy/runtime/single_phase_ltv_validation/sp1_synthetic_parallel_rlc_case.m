function data = sp1_synthetic_parallel_rlc_case(config, opts)
%SP1_SYNTHETIC_PARALLEL_RLC_CASE Synthetic stationary parallel R-L-C load.
%
%   i = G v + Gamma primitive(v) + C dv/dt,   Gamma = 1/L.

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
dv = vEval.values{2};
qv = vEval.primitive;
i = truth.G * v + truth.Gamma * qv + truth.C * dv;

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = ['parallel_rlc_' voltageKind];
data.label = ['Single-phase parallel R-L-C | ' voltageKind];
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'qv', qv, 'dv', dv, 'i', i);
data.truth = truth;
data.truth.law = 'parallel_rlc';
data.noiseSnrDb = noiseSnrDb;
data.voltageKind = voltageKind;
data.voltageComponents = voltageComponents;
data.harmonics = unique(voltageComponents(:, 2)).';
data.power = struct('measured', vN(:) .* iN(:), 'clean', v(:) .* i(:));

end

function truth = default_truth()
truth = struct();
truth.R = 80;
truth.G = 1 / truth.R;
truth.L = 0.120;
truth.Gamma = 1 / truth.L;
truth.C = 25e-6;
truth.voltageRms = 120;
end

function truth = complete_truth(truth)
if ~isfield(truth, 'G') && isfield(truth, 'R')
    truth.G = 1 / truth.R;
end
if ~isfield(truth, 'R') && isfield(truth, 'G')
    truth.R = 1 / truth.G;
end
if ~isfield(truth, 'Gamma') && isfield(truth, 'L')
    truth.Gamma = 1 / truth.L;
end
if ~isfield(truth, 'L') && isfield(truth, 'Gamma')
    truth.L = 1 / truth.Gamma;
end
end

function components = default_voltage_components(kind, truth)
vpk = truth.voltageRms * sqrt(2);
switch kind
    case 'sine'
        components = [vpk, 1, 0.10];
    case 'multisine'
        components = [ ...
            vpk,        1,  0.10; ...
            0.10 * vpk, 3, -0.45; ...
            0.06 * vpk, 5,  1.15; ...
            0.04 * vpk, 7, -1.50];
    otherwise
        error('Unknown voltage kind "%s".', kind);
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
