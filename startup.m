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
if isunix %&& strcmp(char(java.lang.System.getProperty('user.name')), 'dmarx')
    try
        % for all hcit and aftac machines:
        PYENV = pyenv('Version', '/usr/local/bin/python3.7', "ExecutionMode", "InProcess")
        %PYENV = pyenv('Version', '/bin/python3.8', "ExecutionMode", "InProcess")
        py_path = py.sys.path;
        py_path.append('/home/dmarx/HCIT/hcim_mkland3/hcim')
        py_path.append('/home/dmarx/HCIT/hcim_mkland3/hcim/extern')
        %ds9 = @py.ly.util.ds9bdk.ds9fits;
    catch ME
        warning('failed to start python');
        disp(ME.getReport);
    end
else
    try
        PYENV = pyenv('Version', 'C:\Users\dmarx\AppData\Local\conda\conda\envs\py310\python.exe', "ExecutionMode", "InProcess")
    catch ME
        warning('failed to start python');
        disp(ME.getReport);
    end
end

% graphics
set(0,'defaultAxesfontsize',14);

more on
