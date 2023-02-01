mytools = dir('*toolbox');
for i = 1:length(mytools),
   tooldir = [pwd '/' mytools(i).name];
   addpath(tooldir);
   %pdep = genpath(tooldir);
   %addpath(tooldir,pdep);
end

addpath(pwd);

clear;

constants;
unitsdefinitions;

% % python path
if isunix,
    try
        % for all hcit and aftac machines:
        PYENV = pyenv('Version', '/usr/local/bin/python3.7')
        py_path = py.sys.path;
        py_path.append('/home/dmarx/HCIT/hcim_mkland3/hcim')
        py_path.append('/home/dmarx/HCIT/hcim_mkland3/hcim/extern')
        %ds9 = @py.ly.util.ds9bdk.ds9fits;
    catch ME
        warning('failed to start python');
        disp(ME.getReport);
    end
end

% graphics
set(0,'defaultAxesfontsize',14);

more on
