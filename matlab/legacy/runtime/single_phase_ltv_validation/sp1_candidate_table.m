function tableOut = sp1_candidate_table(candidates)
%SP1_CANDIDATE_TABLE Convert candidate structs into a ranking table.

if isempty(candidates)
    tableOut = cell2table(cell(0, 11), 'VariableNames', variable_names());
    return;
end

rows = cell(numel(candidates), 11);
for idx = 1:numel(candidates)
    c = candidates(idx);
    rows(idx, :) = {c.Name, c.Family, c.PhysicalMeaning, c.Accepted, ...
        c.Confidence, c.Coverage, c.ResidualOrScore, c.StabilityCV, ...
        c.Complexity, c.Parameters, c.Reason};
end

tableOut = cell2table(rows, 'VariableNames', variable_names());

end

function names = variable_names()
names = {'Candidate', 'Family', 'PhysicalMeaning', 'Accepted', ...
    'Confidence', 'Coverage', 'ResidualOrScore', 'StabilityCV', ...
    'Complexity', 'Parameters', 'Reason'};
end
