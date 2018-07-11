function [Ex, Ey, bd] = readZBF(filename)
% [Ex, Ey, bd] = readZBF(filename)
% Ex = complex field
% Ey = complex field (=0 if unpolarized)
% bd is struct with fields:
%   Nx, Ny = size of Ex and Ey
%   polflag = 0 => unpolarized, Ey = 0, 1 => polarized
%   dx, dy = sample spacing
%   zz = position relative to pilot beam waist
%   zr = Rayleigh distance for the pilot beam
%   ww = beam waist of the pilot beam
%   lam = wavelength [lens units]

U = CConstants;

fid = fopen(filename);
if fid == -1, error(['error opening ' filename]); end

% read first 9 parameter integers
nn = 9; [A, cnt] = fread(fid,nn,'int32'); if cnt ~= nn, error('read error'); end

formatid = A(1); % currently 1
if formatid ~= 1, error('format error'); end

Nx = A(2); Ny = A(3);
polflag = A(4);  % 0 => unpolarized, 1 => polarized

% Units, 0 for mm, 1 for cm, 2 for in, 3 for meters.
switch A(5), 
    case 0,
        units = U.MM;
    case 1,
        units = U.CM;
    case 2,
        units = U.IN;
    case 3,
        units = 1;
    otherwise,
        error(['unknown units flag: ' num2str(A(5))]);
end

% A(6:9) are unused

% read next 11 parameter doubles
nn = 20; [A, cnt] = fread(fid,nn,'double'); if cnt ~= nn, error('read error'); end

bd = struct(...
    'nx', double(Nx) ...
    ,'ny', double(Ny) ...
    ,'polflag', polflag == 1 ... % boolean, true if polarized
    ,'dx', A(1)*units ... % 1 double: The x direction spacing between points.
    ,'dy', A(2)*units ... % 1 double: The y direction spacing between points.
    ,'zwx', A(3)*units ... % 1 double: The z position relative to the pilot beam waist, x direction.
    ,'zrx', A(4)*units ... % 1 double: The Rayleigh distance for the pilot beam, x direction.
    ,'wwx', A(5)*units ... % 1 double: The waist in lens units of the pilot beam, x direction.
    ,'zwy', A(6)*units ... % 1 double: The z position relative to the pilot beam waist, y direction.
    ,'zry', A(7)*units ... % 1 double: The Rayleigh distance for the pilot beam, y direction.
    ,'wwy', A(8)*units ... % 1 double: The waist in lens units of the pilot beam, y direction.
    ,'wavelength', A(9)*units ... % 1 double: The wavelength in lens units of the beam in the current medium.
    ,'ior', A(10) ... % 1 double: The index of refraction in the current medium.
    ,'FiberRecvEff', A(11) ... % 1 double: The receiver efficiency. Zero if fiber coupling is not computed.
    ,'FiberSysEff', A(12) ... % 1 double: The system efficiency. Zero if fiber coupling is not computed.
);
% 8 doubles: Currently unused, may be any value.

[A, cnt] = fread(fid,[2*Nx,Ny],'double'); if cnt ~= 2*Nx*Ny, error('read error'); end
Ex = A(1:2:end,:)' + 1i*A(2:2:end,:)';

if polflag && ~feof(fid),
   [A, cnt] = fread(fid,[2*Nx,Ny],'double'); if cnt ~= 2*Nx*Ny, error('read error'); end
   Ey = A(1:2:end,:)' + 1i*A(2:2:end,:)';
else
   Ey = 0;
end

fclose(fid);

end % main
