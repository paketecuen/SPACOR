function added = add_legacy_paths()
%ADD_LEGACY_PATHS Add legacy validation folders required by migrated APIs.

paths = spacor.core.legacy_paths();
added = {};
for idx = 1:numel(paths)
    p = paths{idx};
    if isfolder(p)
        addpath(p);
        added{end + 1, 1} = p; %#ok<AGROW>
    else
        warning('spacor:missingLegacyPath', ...
            'Legacy validation path not found: %s', p);
    end
end
end

