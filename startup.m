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
    % for all hcit and aftac machines:
    PYENV = pyenv('Version', '/usr/local/anaconda3/bin/python3.7')
end
    
more on

cd /home/hcit/OMC/OMC_MSWC/falco_testbed_run200
startup
