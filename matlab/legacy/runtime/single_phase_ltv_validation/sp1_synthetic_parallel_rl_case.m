function data = sp1_synthetic_parallel_rl_case(config, opts)
%SP1_SYNTHETIC_PARALLEL_RL_CASE Synthetic stationary parallel R-L load.
%
%   i = G v + Gamma primitive(v),   Gamma = 1/L.

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

truth = nv3_get_option(opts, 'truth', default_truth());
voltageKind = char(nv3_get_option(opts, 'voltageKind', 'sine'));
voltageComponents = nv3_get_option(opts, 'voltageComponents', ...
    default_voltage_components(voltageKind, truth));

vEval = nv3_harmonic_eval(voltageComponents, t, f0, 0);
v = vEval.values{1};
qv = vEval.primitive;
i = truth.G * v + truth.Gamma * qv;

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = ['parallel_rl_' voltageKind];
data.label = ['Single-phase parallel R-L | ' voltageKind];
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'qv', qv, 'i', i);
data.truth = truth;
data.truth.law = 'parallel_rl';
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
truth.voltageRms = 120;
end

function components = default_voltage_components(kind, truth)
vpk = truth.voltageRms * sqrt(2);
switch kind
    case 'sine'
        components = [vpk, 1, 0.15];
    case 'multisine'
        components = [ ...
            vpk,        1,  0.15; ...
            0.10 * vpk, 3, -0.40; ...
            0.05 * vpk, 5,  1.10; ...
            0.03 * vpk, 7, -1.30];
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
