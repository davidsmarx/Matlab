mytools = dir('*toolbox');
for i = 1:length(mytools),
   tooldir = [pwd '/' mytools(i).name];
   addpath(tooldir);
   %pdep = genpath(tooldir);
   %addpath(tooldir,pdep);
end

addpath(pwd); % \Tamar\Matlab

clear;

constants;
unitsdefinitions;

% % python path
if isunix,
    try
        % for all hcit and aftac machines:
        PYENV = pyenv('Version', '/usr/local/bin/python3.7')
    catch ME
        warning('failed to start python');
        disp(ME.getReport);
    end
end
    
more on
