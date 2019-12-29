function LS_startup()

file_path = mfilename('fullpath');
tmp = strfind(file_path,'LS_startup');
file_path = file_path(1:(tmp(end)-1));

% % Folder for all utility functions
% addpath([file_path 'util']);

% Folder for all source files recursively
addpath(genpath([file_path 'src']));

% 
% Folder for all source files recursively
addpath(genpath([file_path 'examples']));

end