function candidate = sp1_candidate_result(fields)
%SP1_CANDIDATE_RESULT Normalize a physical characterization candidate.

candidate = struct();
candidate.Name = get_field(fields, 'Name', '');
candidate.Family = get_field(fields, 'Family', '');
candidate.PhysicalMeaning = get_field(fields, 'PhysicalMeaning', '');
candidate.Accepted = logical(get_field(fields, 'Accepted', false));
candidate.Confidence = get_field(fields, 'Confidence', 0);
candidate.Coverage = get_field(fields, 'Coverage', NaN);
candidate.ResidualOrScore = get_field(fields, 'ResidualOrScore', NaN);
candidate.StabilityCV = get_field(fields, 'StabilityCV', NaN);
candidate.Complexity = get_field(fields, 'Complexity', NaN);
candidate.Parameters = get_field(fields, 'Parameters', '');
candidate.Reason = get_field(fields, 'Reason', '');
candidate.Method = get_field(fields, 'Method', []);
candidate.Diagnostics = get_field(fields, 'Diagnostics', struct());

end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
