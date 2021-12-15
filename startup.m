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

PYENV = pyenv('Version', 'C:\ProgramData\Anaconda3\python.exe')
% % python path
% [pyver, pyexe, isloaded] = pyversion;
% if ~isempty(pyver),
%     pypath = py.sys.path;
%     pypath.append('/home/dmarx/src/Falco-jpl/FALCO-python/falco')
%     pypath.append('/home/dmarx/src/Falco-jpl/FALCO-python/falco/models')
%     pypath.append('/home/dmarx/src/python_toolbox')
% else,
%     warning('Python not found');
% end


more on
