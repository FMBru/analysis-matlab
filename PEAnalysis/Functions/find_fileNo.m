function [out] = find_fileNo(path)
%FIND_FILENO Summary of this function goes here
%   Detailed explanation goes here

    dirdata = dir(path);                    % list contents of the directory
    namelast = dirdata(end).name;           % find last file name 
    out=str2num(namelast(end-8:end-4));     % find last index

end

