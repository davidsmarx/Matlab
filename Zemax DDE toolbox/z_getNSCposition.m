function [x, y, z, x_tilt, y_tilt, z_tilt, material] =...
    z_getNSCposition(zchan, nsurf, nobject)
% [x, y, z, x_tilt, y_tilt, z_tilt] = z_getNSCposition(zchan, nsurf, nobject)
% sPos = z_getNSCposition(zchan, nsurf, nobject)
%
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC

UU = CConstants;

cmdstr = sprintf('GetNSCPosition,%d,%d',nsurf,nobject);
opvalstr = ddereq(zchan,cmdstr,[1 1]);

opvals = textscan(opvalstr,'%f,%f,%f,%f,%f,%f,%s');

x = opvals{1}*UU.MM;
y = opvals{2}*UU.MM;
z = opvals{3}*UU.MM;
x_tilt = opvals{4}*UU.P;
y_tilt = opvals{5}*UU.P;
z_tilt = opvals{6}*UU.P;
material = opvals{7};

if nargout == 1,
    sPos = struct('x',x,'y',y,'z',z,...
        'x_tilt',x_tilt,'y_tilt',y_tilt,'z_tilt',z_tilt,...
        'material',material);
    x = sPos;
end

return
