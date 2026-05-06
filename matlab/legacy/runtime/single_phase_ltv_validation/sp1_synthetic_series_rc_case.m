function data = sp1_synthetic_series_rc_case(config, opts)
%SP1_SYNTHETIC_SERIES_RC_CASE Synthetic stationary series R-C load.
%
%   v = R i + Gamma primitive(i),   Gamma = 1/C.

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
currentKind = char(nv3_get_option(opts, 'currentKind', 'sine'));
currentComponents = nv3_get_option(opts, 'currentComponents', ...
    default_current_components(currentKind, truth));

iEval = nv3_harmonic_eval(currentComponents, t, f0, 0);
i = iEval.values{1};
qi = iEval.primitive;
v = truth.R * i + truth.Gamma * qi;

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = ['series_rc_' currentKind];
data.label = ['Single-phase series R-C | ' currentKind];
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'i', i, 'qi', qi);
data.truth = truth;
data.truth.law = 'series_rc';
data.noiseSnrDb = noiseSnrDb;
data.currentKind = currentKind;
data.currentComponents = currentComponents;
data.harmonics = unique(currentComponents(:, 2)).';
data.power = struct('measured', vN(:) .* iN(:), 'clean', v(:) .* i(:));

end

function truth = default_truth()
truth = struct();
truth.R = 12;
truth.C = 80e-6;
truth.Gamma = 1 / truth.C;
truth.currentRms = 8;
end

function components = default_current_components(kind, truth)
ipk = truth.currentRms * sqrt(2);
switch kind
    case 'sine'
        components = [ipk, 1, -0.20];
    case 'multisine'
        components = [ ...
            ipk,        1, -0.20; ...
            0.10 * ipk, 3,  0.55; ...
            0.05 * ipk, 5, -1.00; ...
            0.03 * ipk, 7,  1.60];
    otherwise
        error('Unknown current kind "%s".', kind);
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
