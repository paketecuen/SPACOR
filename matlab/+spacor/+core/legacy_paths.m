function paths = legacy_paths()
%LEGACY_PATHS Return validated runtime folders bundled with this release.
%
% Public APIs live under +spacor.  Some low-level routines from the validated
% research runtime are kept under matlab/legacy/runtime until they are fully
% refactored into package namespaces.  No private data roots are used here.

labRoot = spacor.core.lab_root();

singlePhaseRuntime = fullfile(labRoot, 'legacy', 'runtime', ...
    'single_phase_ltv_validation');

threeWireRuntime = fullfile(labRoot, 'legacy', 'runtime', ...
    'three_wire_validation');

paths = { ...
    singlePhaseRuntime, ...
    threeWireRuntime ...
    };
end
