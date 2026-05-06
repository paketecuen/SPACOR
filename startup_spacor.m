function root = startup_spacor()
%STARTUP_SPACOR Convenience startup from the repository root.
%
% Usage from MATLAB:
%   cd /path/to/SPACOR
%   startup_spacor

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'matlab'));
root = startup_spacor_lab();
end
