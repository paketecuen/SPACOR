function candidate = sp1_candidate_lti_series_rl(data, config, opts, decisionOpts)
%SP1_CANDIDATE_LTI_SERIES_RL Candidate wrapper for series R-L.

if nargin < 3
    opts = struct();
end
if nargin < 4
    decisionOpts = struct();
end

base = sp1_candidate_windowed_rl(data, config, opts, decisionOpts);
base.Name = 'lti_series_rl';
base.Family = 'lti_windowed';
base.PhysicalMeaning = 'series R-L: v = R i + L di/dt';
base.Complexity = 2;
candidate = base;

end
