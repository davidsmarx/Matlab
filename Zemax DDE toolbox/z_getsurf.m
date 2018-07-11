function surfdata = z_getsurf(zchan,surf)
% surfdata = z_getsurf(zchan,nsurf)
%
% nsurf = integer surface #, or surf struct
%
% surfdata is a surf struct with the following fields:
% nsurf, type, curv, thick, glass, radius (aperture), conic, parms, coat

global MM UM;

% parse input
switch class(surf),
    case 'double',
        nsurf = surf;
    case 'struct',
        nsurf = surf.nsurf;
    otherwise,
        error('unknown surf type');
end

surfstr = ddereq(zchan,['GetSurface,' num2str(nsurf)],[1 1]);
surfstr = strtrim(surfstr); % removes leading and trailing white space

seps = findstr(surfstr,',');

surfdata.nsurf = nsurf;
surfdata.type = surfstr(1:seps(1)-1);
surfdata.curv = str2num(surfstr(seps(1)+1:seps(2)-1))*(1/MM);
surfdata.thick = str2num(surfstr(seps(2)+1:seps(3)-1))*MM;
surfdata.glass = surfstr(seps(3)+1:seps(4)-1);
surfdata.radius = str2num(surfstr(seps(4)+1:seps(5)-1))*MM;
surfdata.conic = str2num(surfstr(seps(5)+1:seps(6)-1));
surfdata.parms = sscanf(surfstr(seps(6)+1:end),'%e,',8);
surfdata.coat = surfstr(seps(end)+1:end);

return