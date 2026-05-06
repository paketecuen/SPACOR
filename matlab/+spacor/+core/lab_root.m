function root = lab_root()
%LAB_ROOT Return the canonical SPACOR MATLAB package root folder.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
