function report = sp1_characterize_physical_load(data, config, opts)
%SP1_CHARACTERIZE_PHYSICAL_LOAD Physical model selection for one-phase loads.
%
% This is the high-level API intended for black-box characterization.  It
% runs candidate physical models and returns an engineering-oriented report:
%
%   1. windowed LTI RLC candidates;
%   2. switching/topology segmentation;
%   3. conservative nonlinear storage lambda=f(i);
%   4. hysteretic B-H loop.
%
% The goal is not to estimate parameters for their own sake, but to decide
% which physical representation is useful for monitoring/diagnosis.

if nargin < 2 || isempty(config)
    config = sp1_default_config();
end
if nargin < 3
    opts = struct();
end

data = normalize_data(data, config);
opts = with_defaults(opts, config);
candidateMode = lower(string(opts.candidateMode));

ltiCandidateResults = [
    sp1_candidate_lti_series_rl(data, config, opts.seriesRL, opts)
    sp1_candidate_lti_parallel_rl(data, config, opts.parallelRL, opts)
    sp1_candidate_lti_parallel_series_rl(data, config, ...
    opts.parallelSeriesRL, opts)
    sp1_candidate_lti_series_rc(data, config, opts.seriesRC, opts)
    sp1_candidate_lti_parallel_rc(data, config, opts.parallelRC, opts)
    sp1_candidate_lti_capacitor(data, config, opts.capacitor, opts)
    sp1_candidate_lti_series_rlc(data, config, opts.seriesRLC, opts)
    sp1_candidate_lti_parallel_rlc(data, config, opts.parallelRLC, opts)];

switch candidateMode
    case "lti"
        candidateResults = ltiCandidateResults;
    case "full"
        candidateResults = [
            ltiCandidateResults
            sp1_candidate_switching_affine_rl(data, config, opts.switching)
            sp1_candidate_flux_curve(data, config, opts.fluxCurve, opts)
            sp1_candidate_bh_loop(data, opts.bhLoop, opts)];
    otherwise
        error('Unknown candidateMode "%s". Use "lti" or "full".', ...
            candidateMode);
end

candidates = sp1_candidate_table(candidateResults);
recommendation = sp1_select_candidate_recommendation(candidates, opts);

methods = struct();
methods.seriesRL = ltiCandidateResults(1).Method;
methods.windowedRL = ltiCandidateResults(1).Method;
methods.parallelRL = ltiCandidateResults(2).Method;
methods.parallelSeriesRL = ltiCandidateResults(3).Method;
methods.seriesRC = ltiCandidateResults(4).Method;
methods.parallelRC = ltiCandidateResults(5).Method;
methods.capacitor = ltiCandidateResults(6).Method;
methods.seriesRLC = ltiCandidateResults(7).Method;
methods.parallelRLC = ltiCandidateResults(8).Method;
methods.switching = [];
methods.fluxCurve = [];
methods.bhLoop = [];
if candidateMode == "full"
    methods.switching = candidateResults(9).Method;
    methods.fluxCurve = candidateResults(10).Method;
    methods.bhLoop = candidateResults(11).Method;
end

report = struct();
report.input = input_summary(data);
report.candidates = candidates;
report.candidateResults = candidateResults;
report.recommendation = recommendation;
report.methods = methods;
report.flow = characterization_flow();
report.options = opts;

end

function data = normalize_data(data, config)
required = {'v', 'i'};
for idx = 1:numel(required)
    if ~isfield(data, required{idx})
        error('Data must contain field "%s".', required{idx});
    end
end
data.v = data.v(:);
data.i = data.i(:);
if ~isfield(data, 'fs') || isempty(data.fs)
    data.fs = config.fs;
end
if ~isfield(data, 'f0') || isempty(data.f0)
    data.f0 = config.f0;
end
if ~isfield(data, 't') || isempty(data.t)
    data.t = (0:numel(data.v) - 1).' / data.fs;
else
    data.t = data.t(:);
end
if numel(data.i) ~= numel(data.v) || numel(data.t) ~= numel(data.v)
    error('Data fields v, i, and t must have the same length.');
end
if ~isfield(data, 'label')
    data.label = 'single-phase black-box load';
end
end

function opts = with_defaults(opts, config)
opts.candidateMode = nv3_get_option(opts, 'candidateMode', 'full');
opts.residualGate = nv3_get_option(opts, 'residualGate', 2e-2);
opts.coverageGate = nv3_get_option(opts, 'coverageGate', 0.80);
opts.stabilityGate = nv3_get_option(opts, 'stabilityGate', 0.15);
opts.fluxScoreGate = nv3_get_option(opts, 'fluxScoreGate', 5e-2);
opts.bhResidualGate = nv3_get_option(opts, 'bhResidualGate', 2e-2);
opts.bhImprovementFactor = nv3_get_option(opts, 'bhImprovementFactor', 0.50);
opts.ltiEquivalenceFactor = nv3_get_option(opts, 'ltiEquivalenceFactor', 1.5);
opts.ltiEquivalenceAbsTol = nv3_get_option(opts, 'ltiEquivalenceAbsTol', 2e-3);
if ~isfield(opts, 'seriesRL') || isempty(opts.seriesRL)
    if isfield(opts, 'windowedRL') && ~isempty(opts.windowedRL)
        opts.seriesRL = opts.windowedRL;
    else
        opts.seriesRL = struct();
    end
end
opts.seriesRL.windowCycles = nv3_get_option(opts.seriesRL, ...
    'windowCycles', config.windowCycles);
opts.seriesRL.hopCycles = nv3_get_option(opts.seriesRL, ...
    'hopCycles', config.hopCycles);
opts.seriesRL.harmonics = nv3_get_option(opts.seriesRL, ...
    'harmonics', config.harmonics);
opts.windowedRL = opts.seriesRL;
if ~isfield(opts, 'parallelRL') || isempty(opts.parallelRL)
    opts.parallelRL = struct();
end
opts.parallelRL.windowCycles = nv3_get_option(opts.parallelRL, ...
    'windowCycles', max(config.windowCycles, 0.5));
opts.parallelRL.hopCycles = nv3_get_option(opts.parallelRL, ...
    'hopCycles', config.hopCycles);
opts.parallelRL.harmonics = nv3_get_option(opts.parallelRL, ...
    'harmonics', config.harmonics);
if ~isfield(opts, 'parallelSeriesRL') || isempty(opts.parallelSeriesRL)
    opts.parallelSeriesRL = struct();
end
opts.parallelSeriesRL.windowCycles = nv3_get_option(opts.parallelSeriesRL, ...
    'windowCycles', max(config.windowCycles, 0.5));
opts.parallelSeriesRL.hopCycles = nv3_get_option(opts.parallelSeriesRL, ...
    'hopCycles', config.hopCycles);
opts.parallelSeriesRL.harmonics = nv3_get_option(opts.parallelSeriesRL, ...
    'harmonics', config.harmonics);
if ~isfield(opts, 'seriesRC') || isempty(opts.seriesRC)
    opts.seriesRC = struct();
end
opts.seriesRC.windowCycles = nv3_get_option(opts.seriesRC, ...
    'windowCycles', max(config.windowCycles, 0.5));
opts.seriesRC.hopCycles = nv3_get_option(opts.seriesRC, ...
    'hopCycles', config.hopCycles);
opts.seriesRC.harmonics = nv3_get_option(opts.seriesRC, ...
    'harmonics', config.harmonics);
if ~isfield(opts, 'parallelRC') || isempty(opts.parallelRC)
    opts.parallelRC = struct();
end
opts.parallelRC.windowCycles = nv3_get_option(opts.parallelRC, ...
    'windowCycles', config.windowCycles);
opts.parallelRC.hopCycles = nv3_get_option(opts.parallelRC, ...
    'hopCycles', config.hopCycles);
opts.parallelRC.harmonics = nv3_get_option(opts.parallelRC, ...
    'harmonics', config.harmonics);
if ~isfield(opts, 'capacitor') || isempty(opts.capacitor)
    opts.capacitor = struct();
end
opts.capacitor.windowCycles = nv3_get_option(opts.capacitor, ...
    'windowCycles', config.windowCycles);
opts.capacitor.hopCycles = nv3_get_option(opts.capacitor, ...
    'hopCycles', config.hopCycles);
opts.capacitor.harmonics = nv3_get_option(opts.capacitor, ...
    'harmonics', config.harmonics);
if ~isfield(opts, 'seriesRLC') || isempty(opts.seriesRLC)
    opts.seriesRLC = struct();
end
opts.seriesRLC.windowCycles = nv3_get_option(opts.seriesRLC, ...
    'windowCycles', max(config.windowCycles, 0.5));
opts.seriesRLC.hopCycles = nv3_get_option(opts.seriesRLC, ...
    'hopCycles', config.hopCycles);
opts.seriesRLC.harmonics = nv3_get_option(opts.seriesRLC, ...
    'harmonics', config.harmonics);
if ~isfield(opts, 'parallelRLC') || isempty(opts.parallelRLC)
    opts.parallelRLC = struct();
end
opts.parallelRLC.windowCycles = nv3_get_option(opts.parallelRLC, ...
    'windowCycles', max(config.windowCycles, 0.5));
opts.parallelRLC.hopCycles = nv3_get_option(opts.parallelRLC, ...
    'hopCycles', config.hopCycles);
opts.parallelRLC.harmonics = nv3_get_option(opts.parallelRLC, ...
    'harmonics', config.harmonics);
if ~isfield(opts, 'fluxCurve') || isempty(opts.fluxCurve)
    opts.fluxCurve = struct();
end
opts.fluxCurve.degree = nv3_get_option(opts.fluxCurve, 'degree', 9);
opts.fluxCurve.rGrid = nv3_get_option(opts.fluxCurve, ...
    'rGrid', linspace(0, 2, 401));
if ~isfield(opts, 'bhLoop') || isempty(opts.bhLoop)
    opts.bhLoop = struct();
end
opts.bhLoop.rGrid = nv3_get_option(opts.bhLoop, ...
    'rGrid', linspace(0, 2, 401));
opts.bhLoop.includeHysteresis = nv3_get_option(opts.bhLoop, ...
    'includeHysteresis', true);
if ~isfield(opts, 'switching') || isempty(opts.switching)
    opts.switching = struct();
end
end

function summary = input_summary(data)
summary = struct('label', data.label, 'samples', numel(data.t), ...
    'fs', data.fs, 'f0', data.f0, ...
    'duration', data.t(end) - data.t(1), ...
    'vRms', sqrt(mean(data.v .^ 2)), ...
    'iRms', sqrt(mean(data.i .^ 2)));
end

function flow = characterization_flow()
flow = { ...
    'Try simple LTI windowed energetic equivalents.'; ...
    'Prefer switching if ON/OFF topology is observable.'; ...
    'If LTI is unstable or residual is high, try lambda=f(i).'; ...
    'If lambda=f(i) is not single-valued, try a hysteretic B-H loop.'; ...
    'If no model passes, report missing physics or insufficient excitation.'};
end
