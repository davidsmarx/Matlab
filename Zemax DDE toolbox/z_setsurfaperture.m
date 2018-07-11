function status = z_setsurfaperture(zchan, surf, strAperture)
% status = z_setsurfaperture(zchan, surf, strAperture)
% 
% surf is a surface number or a surface struct as returned by z_getsurf
%
% strAperture is a struct with the following fields:
%   type
%   min
%   max
%   decenter = [x, y]
%   aperturefilename
%
% strAperture = z_getsurfaperture with no input arguments returns a default
% struct

% set units
global MM;

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

% SetAperture, surf, type, min, max, xdecenter, ydecenter, aperturefile
cmdstr = ['SetAperture, ' ...
    num2str(nsurf) ', ' ...
    num2str(strAperture.type) ', ' ...
    num2str(strAperture.min/MM) ', ' ...
    num2str(strAperture.max/MM) ', ' ...
    num2str(strAperture.decenter(1)/MM) ', ' ... 
    num2str(strAperture.decenter(2)/MM) ', ' ...
    '"' strAperture.aperturefilename '"' ...
    ];

status = ddereq(zchan, cmdstr, [1 1]);
    

return