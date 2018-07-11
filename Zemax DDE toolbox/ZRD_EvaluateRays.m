function S = ZRD_EvaluateRays( rays, sOptions )
% S = ZRD_EvaluateRays( rays, sOptions )
%
% rays is an array of struct (returned by mexReadZRD and ZRD_ReadZRD):
%
% 	const char *field_names[] = { "status"
% 	,"level"
% 	,"hit_object"
% 	,"hit_face"
% 	,"unused"
% 	,"in_object"
% 	,"parent"
% 	,"storage"
% 	,"xybin", "lmbin"
% 	,"index", "starting_phase"
% 	,"x", "y", "z"
% 	,"l", "m", "n"
% 	,"nx", "ny", "nz"
% 	,"path_to", "intensity"
% 	,"phase_of", "phase_at"
% 	,"exr", "exi", "eyr", "eyi", "ezr", "ezi"
%   ,"wavelength"};
%
% 
% sOptions = struct(... (use sParms returned by z_GetDetectorViewer)
%     'DetectorSize', [0.009200000000000 0.009200000000000] ...
%     ,'DetectorPixels', [64 64] ...
%     );
%
% IrrMap(Nx,Ny) is [W/m^2]

sOptions = ValidateOptions(sOptions);

% pixel area in cm^2:
%pixel_area_cm2 = (9.2*MM/CM/64).^2;

Nx = sOptions.DetectorPixels(1);
Ny = sOptions.DetectorPixels(2);
dx = sOptions.DetectorSize(1)/Nx; % 9.2*MM/N;
dy = sOptions.DetectorSize(2)/Ny;
x = dx*(-(Nx/2-1/2):(Nx/2-1/2))';
y = dy*(-(Ny/2-1/2):(Ny/2-1/2))';

% the zemax rays are x-major
[Y, X] = meshgrid(y, x);
% from the zemax manual: For the Detector Rectangle, the pixel numbers start at 1 in the lower left (-x, -y) corner of the rectangle. The
% pixel numbers increase along the +x direction first, and then move up one row in the +y direction, and start over
% at the -x edge. For a detector with nx by ny pixels, the pixel numbers are 1 through nx on the bottom row, and
% then nx+1 through 2nx on the next row above the bottom, until the last pixel (number nx*ny) is reached in the
% upper right (+x, +y) corner. Another way of stating the pixel numbering is that the x index changes the fastest,
% and then the y index. check: 
% % [X(1) Y(1)] is lower left:
% % [X(N) Y(N)] is lower right:
% % [X(N+1) Y(N+1)] is left side, up one row (Y)
% % [X(N*N) Y(N*N)] is upper right:

% calculate phase at each of the four pixels
% the detector object is rotated 180 about y, but rays.x, rays.y are global
% coords, but detbin is local. So need to un-rotate local X(bin) back to
% global coords.

% calculate incoherent irradiance
numrays = length(rays);
[IrrMap, AmpMap, FldMap] = deal(zeros(size(X)));
for ii = 1:numrays,
    
    thisray = rays(ii);
        
    % looks like xybin is zero-offset, not one-offset
    detbin = thisray.xybin+1;
    xp = X(detbin); yp = Y(detbin);
    % the detector object is rotated 180 about y, rays.x, rays.y are global
    % coords, but detbin is local. So need to rotate global rays.x,y about
    % y:
    xr = -thisray.x;        yr = thisray.y;

    % double check that ray x,y is closest to bin center given by xybin
    if abs(xr-xp) > dx/2 || abs(yr-yp) > dy/2,
        error(['x,y outside bin, ray = ' num2str(ii)]);
    end

    % distance from pixel center
    ddx = xr-xp;
    ddy = yr-yp;
    
    dirx = sign(ddx);
    diry = sign(ddy);

    % the other three closest detector pixels
    detbin_dx = detbin + dirx;
    detbin_dy = detbin + Nx*diry;
    detbin_dxdy = detbin + dirx + Nx*diry;
    
    DetBin = [
        detbin detbin_dx
        detbin_dy detbin_dxdy
        ];
    
    % calculate phase at each of the four pixels
    k = 2*pi./thisray.wavelength;
        
    % ray.phase_at from the zrd file is already calculated at the pixel
    % center by zemax. see e-mail exchange with zemax support
    Pha = thisray.phase_at + ... 
        k.* ( (-X(detbin) - -X(DetBin)).*thisray.l + (Y(detbin) - Y(DetBin)).*thisray.m );

    % calculate field amplitude (=sqrt(intensity)) at each of the four
    % pixels

    % bi-linear interpolation to apportion ray value to four pixels
    ax = abs(ddx)./dx;
    ay = abs(ddy)./dy;

    
    Amp = sqrt([
        (1-ax)*(1-ay) ax*(1-ay)
        (1-ax)*ay     ax*ay
        ] .* thisray.intensity );
    
    Field = Amp .* exp(1i*Pha);

    % sum the quantities for these pixels
    IrrMap(DetBin) = IrrMap(DetBin) + Amp.*Amp;
    AmpMap(DetBin) = AmpMap(DetBin) + Amp;
    FldMap(DetBin) = FldMap(DetBin) + Field;
    
end % for each ray

% calculate according to equations in zemax manual
CohIrr = IrrMap .* abs(FldMap).^2 ./ AmpMap.^2;
CohPha = angle(FldMap);

% from zemax row-major to matlab column-major
% Watts to Watts/m^2, DetectorViewer output is W/cm^2
S.IrrMap = (IrrMap.')./( dx*dy );
S.CohIrr = (CohIrr.')./( dx*dy );
S.CohPha = (CohPha.');
S.Field  = (FldMap.')./sqrt( dx*dy );
S.x = x;
S.y = y;
S.X = (X.');
S.Y = (Y.');

% figure, imagesc(IrrMap), colormap(flipud(gray)), colorb ar
% axis image

end % main

function sOptions = ValidateOptions(sIn)

% default:
sOptions = struct(...
    'DetectorSize', [0.009200000000000 0.009200000000000] ...
    ,'DetectorPixels', [64 64] ...
    );

fnames = fieldnames(sIn);
for ii = 1:length(fnames),
    sOptions.(fnames{ii}) = sIn.(fnames{ii});
end

end % ValidateOptions

