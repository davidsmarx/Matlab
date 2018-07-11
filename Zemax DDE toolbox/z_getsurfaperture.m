function strAperture = z_getsurfaperture(zchan,surf)
% strAperture = z_getsurfaperture(zchan,surf)
% 
% surf is a surface number or a surface struct as returned by z_getsurf
%
% strAperture is a struct with the following fields:
%   type
%   min radius
%   max radius
%   decenter = [x, y]
%   aperturefilename
%
% strAperture = z_getsurfaperture with no input arguments returns a default
% struct

global MM P;

% initialize return struct
strAperture = struct('type',0,'min',0,'max',0,...
    'decenter',[0 0],'aperturefilename','');

if nargin == 0, return, end

switch class(surf),
    case 'double',
        nsurf = surf;
        strsurf = struct;
    case 'struct',
        nsurf = surf.nsurf;
        strsurf = surf;
    otherwise,
        error('invalid surf input');
end

cmdstr = ['GetAperture,' num2str(nsurf)];
retstr = ddereq(zchan, cmdstr,[1 1]);

CC = textscan(retstr,'%f %f %f %f %f %s','delimiter',',');
[strAperture.type, strAperture.min, strAperture.max,...
    strAperture.decenter(1), strAperture.decenter(2),...
    cellfilename] = deal(CC{:});

strAperture.aperturefilename = cellfilename{1};

% units
strAperture.min = strAperture.min*MM;
strAperture.max = strAperture.max*MM;
strAperture.decenter = strAperture.decenter*MM;

return