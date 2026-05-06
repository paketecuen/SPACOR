function [data, config, truth] = singlephase_series_rl(config)
%SINGLEPHASE_SERIES_RL Deterministic single-phase RL fixture.

if nargin < 1 || isempty(config)
    config = spacor.singlephase.default_config();
end

config.fs = spacor.core.get_option(config, 'fs', 10000);
config.f0 = spacor.core.get_option(config, 'f0', 50);
config.duration = spacor.core.get_option(config, 'duration', 0.20);
config.windowCycles = spacor.core.get_option(config, 'windowCycles', 0.5);
config.hopCycles = spacor.core.get_option(config, 'hopCycles', 0.10);
config.harmonics = spacor.core.get_option(config, 'harmonics', [1 3 5 7]);

t = (0:round(config.duration * config.fs) - 1).' / config.fs;
w0 = 2 * pi * config.f0;
truth = struct('R', 10, 'L', 50e-3, 'model', 'lti_series_rl');
harmonics = [1 3 5 7];
amps = [120 9 5 3];
phases = [0 0.4 -0.7 1.1];

v = zeros(size(t));
i = zeros(size(t));
for idx = 1:numel(harmonics)
    h = harmonics(idx);
    omega = h * w0;
    z = truth.R + 1i * omega * truth.L;
    v = v + amps(idx) * sin(omega * t + phases(idx));
    i = i + amps(idx) / abs(z) * sin(omega * t + phases(idx) - angle(z));
end

data = struct('t', t, 'v', v, 'i', i, 'fs', config.fs, ...
    'f0', config.f0, 'truth', truth, ...
    'label', 'fixture single-phase series RL');
end
