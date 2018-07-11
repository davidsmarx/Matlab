mytools = dir('*toolbox');
for i = 1:length(mytools),
   tooldir = [pwd '/' mytools(i).name];
   addpath(tooldir);
   %pdep = genpath(tooldir);
   %addpath(tooldir,pdep);
end

addpath(pwd); % \Tamar\Matlab

clear;

more on

constants;
unitsdefinitions;

% % hardware
% global IMSMOTION;
% IMSMOTION = 1;
% 
% %Aq_Startup.m   Acqiris Startup file
% % 
% % This file is generated automatically during the installation
% % DO NOT EDIT - File will be overwritten if re-installing or updating
%  
% AqRoot=getenv('AcqirisDxRoot');
% if ~isempty(AqRoot)
%     addpath([AqRoot,'\bin']);
%     addpath([AqRoot,'\MATLAB\mex']); %,[AqRoot,'\MATLAB\mex\help']);
%     if exist([AqRoot,'\MATLAB\mex\functions'],'dir'),
%         addpath([AqRoot,'\MATLAB\mex\functions']);
%     end
% end