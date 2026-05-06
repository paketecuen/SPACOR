function report = sp1_characterize_switching_load(data, config, opts)
%SP1_CHARACTERIZE_SWITCHING_LOAD Switching-aware physical characterization.

if nargin < 2 || isempty(config)
    config = sp1_default_config();
end
if nargin < 3
    opts = struct();
end

if ~isfield(opts, 'segmentation')
    opts.segmentation = struct();
end
if ~isfield(opts, 'identification')
    opts.identification = struct();
end
if ~isfield(opts, 'regime')
    opts.regime = struct();
end

segmentation = sp1_segment_switching_states(data, config, opts.segmentation);
regime = sp1_switching_regime_report(segmentation, data, config, opts.regime);
model = sp1_identify_switched_lti_library(data, segmentation, opts.identification);
accepted = model.summary.Accepted;
confidence = model.summary.Confidence;
reason = model.summary.Reason;

report = struct();
report.model = model.summary.Model;
report.physicalMeaning = model.summary.PhysicalMeaning;
report.accepted = accepted;
report.confidence = confidence;
report.reason = reason;
report.parameters = parameter_string(model.summary);
report.segmentation = segmentation;
report.regime = regime;
report.identification = model;
report.summary = model.summary;

end

function text = parameter_string(summary)
parts = {sprintf('V0=%.6g V', summary.MedianV0)};
if isfinite(summary.MedianR)
    parts{end + 1} = sprintf('R=%.6g ohm', summary.MedianR);
end
if isfinite(summary.MedianL)
    parts{end + 1} = sprintf('L=%.6g H', summary.MedianL);
end
if isfinite(summary.MedianGamma)
    parts{end + 1} = sprintf('Gamma=%.6g', summary.MedianGamma);
end
if isfinite(summary.MedianC)
    parts{end + 1} = sprintf('C=%.6g F', summary.MedianC);
end
text = strjoin(parts, ', ');
end
