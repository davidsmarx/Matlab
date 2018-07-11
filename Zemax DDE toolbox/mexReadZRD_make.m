% make script to compile mexReadZRD

% get this path
mexdir = mfilename('fullpath');
% strip the mfilename from the path to get only the path
mexdir = mexdir(1:end-length(mfilename));

curdir = pwd;
try,
cd(mexdir);
mex mexReadZRD/mexReadZRD/mexReadZRD.cpp mexReadZRD/mexReadZRD/ReadBufferCalcFields.cpp mexReadZRD/mexReadZRD/ProcessRays.cpp;
catch,
    disp('mex failed!');
end
cd(curdir);

