function root = startup_spacor_lab()
%STARTUP_SPACOR_LAB Add the public SPACOR MATLAB package to the MATLAB path.

root = fileparts(mfilename('fullpath'));
repoRoot = fileparts(root);
addpath(root);
addpath(fullfile(repoRoot, 'scripts'));
addpath(fullfile(repoRoot, 'tests'));

fprintf('SPACOR public MATLAB package ready: %s\n', root);
fprintf('APIs are available under +spacor.\n');
end
