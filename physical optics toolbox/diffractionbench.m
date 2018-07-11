function [elapsedtime, timeperthread] = diffractionbench(NT)
% elapsedtime = diffractionbench(NT)
%
% NT = number of threads to use
% NT = 0 => single-threaded interface
% NT = 1 => multi-threaded interface using only one thread

%global MM UM NM;
unitsdefinitions;

if NT == 0,
    TID = 0;
    BID = 0;
    NT = 1;
elseif NT > 0,
    TID = 1:NT;
    BID = 1:NT;
else,
    error('NT must be >= 0');
end

source = struct('nx',1024,'ny',1024,'Gx',15*MM,'Gy',15*MM,'wavelength',1.319*UM);
mask = [3 2 3]*MM;
ccube = IntMetCom_CornerCube;
ccube.gapwidth = [0.5 0 0]*MM;
ccube.rotmatrix = corner_cube_orientations('ew',0,0);

hD = IntMetCom_init;

tic,

for ii = 1:NT,
    hD.ClearThread(TID(ii));
    
    % create gauss beam
    IntMetCom_CreateGausSource(hD,10*MM,source,TID(ii),BID(ii));
    
    for jj = 1:3,
        % propagate Fresnel
        z = source.Gx.^2 ./ (source.nx*source.wavelength) - 1*MM;
        IntMetCom_Propagate(hD,z,TID(ii),BID(ii));
        
        % propagate Fraunhoffer using CZT
        z = source.Gx.^2 ./ (source.nx*source.wavelength) + 1*MM;
        IntMetCom_Propagate(hD,z,TID(ii),BID(ii));
        
        % apply two hole mask
        IntMetCom_ApplyMaskRotate(hD,mask,0,0,0,'hh',TID(ii),BID(ii));
        
        % apply corner cube
        IntMetCom_CornerCube(hD,ccube,TID(ii),BID(ii));
    end
end

for ii = 1:NT,
    hD.ExecuteThread(TID(ii));
end
IntMetCom_WaitForThreads(hD,TID);

elapsedtime = toc;
timeperthread = elapsedtime/NT;

hD.delete;
