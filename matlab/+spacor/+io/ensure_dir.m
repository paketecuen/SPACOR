function pathOut = ensure_dir(pathIn)
%ENSURE_DIR Create a directory if it does not exist.

pathOut = char(pathIn);
if ~exist(pathOut, 'dir')
    mkdir(pathOut);
end
end

