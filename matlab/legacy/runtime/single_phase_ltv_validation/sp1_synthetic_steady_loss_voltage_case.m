function data = sp1_synthetic_steady_loss_voltage_case(config, opts)
%SP1_SYNTHETIC_STEADY_LOSS_VOLTAGE_CASE Voltage-driven steady iron loss RL.
%
% This case emulates iron losses as a stationary dissipative equivalent:
%
%   v = (R + Rfe) i + L di/dt.
%
% Algebraically Rfe occupies the same slot as the dL/dt term in
%
%   v = R i + L di/dt + dL/dt i,
%
% but L is kept periodic/stationary.  This is the correct test when we want
% constant parameters under a stationary sinusoidal or distorted voltage.

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
    default_voltage_components(voltageKind));

vEval = nv3_harmonic_eval(voltageComponents, t, f0, 0);
currentComponents = steady_current_components(voltageComponents, ...
    truth.R + truth.Rfe, truth.L, f0);
iEval = nv3_harmonic_eval(currentComponents, t, f0, 1);
v = vEval.values{1};
i = iEval.values{1};
di = iEval.values{2};

noiseSnrDb = nv3_get_option(opts, 'noiseSnrDb', config.noiseSnrDb);
seed = nv3_get_option(opts, 'randomSeed', config.randomSeed);
[vN, iN] = add_noise(seed, noiseSnrDb, v, i);

data = struct();
data.model = ['series_rl_steady_loss_' voltageKind];
data.label = ['Single-phase steady R + Rfe + L | ' voltageKind];
data.fs = fs;
data.f0 = f0;
data.t = t;
data.v = vN;
data.i = iN;
data.clean = struct('v', v, 'i', i, 'di', di, ...
    'L', truth.L * ones(size(t)), ...
    'Lcoef', truth.L * ones(size(t)), ...
    'Linc', truth.L * ones(size(t)), ...
    'dLdt', truth.Rfe * ones(size(t)), ...
    'lambda', truth.L * i, ...
    'Rloss', truth.Rfe * ones(size(t)), ...
    'Reff', (truth.R + truth.Rfe) * ones(size(t)));
data.truth = truth;
data.truth.Reff = truth.R + truth.Rfe;
data.truth.law = 'steady_loss_equivalent';
data.noiseSnrDb = noiseSnrDb;
data.voltageKind = voltageKind;
data.voltageComponents = voltageComponents;
data.currentComponents = currentComponents;
data.harmonics = unique(voltageComponents(:, 2)).';
data.power = struct('measured', vN(:) .* iN(:), 'clean', v(:) .* i(:));

end

function truth = default_truth()
truth = struct();
truth.R = 0.65;
truth.Rfe = 0.12;
truth.L = 0.080;
end

function components = default_voltage_components(kind)
vpk = 230 * sqrt(2);
switch kind
    case 'sine'
        components = [vpk, 1, 0.10];
    case 'multisine'
        components = [ ...
            vpk,        1,  0.10; ...
            0.08 * vpk, 3, -0.45; ...
            0.05 * vpk, 5,  1.10; ...
            0.03 * vpk, 7, -1.70];
    otherwise
        error('Unknown voltage kind "%s".', kind);
end
end

function currentComponents = steady_current_components(vComponents, r, L, f0)
omega0 = 2 * pi * f0;
currentComponents = zeros(size(vComponents));
for row = 1:size(vComponents, 1)
    amplitude = vComponents(row, 1);
    harmonic = vComponents(row, 2);
    phase = vComponents(row, 3);
    z = r + 1i * harmonic * omega0 * L;
    currentComponents(row, :) = [amplitude / abs(z), harmonic, ...
        phase - angle(z)];
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
