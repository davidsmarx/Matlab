function strsurf = z_getsurftiltdec(zchan,surf)
% strsurf = z_getsurftiltdec(zchan,surf)
% 
% surf is a surface number or a surface struct as returned by z_getsurf
%
% strsurf is a struct with the following fields:
% beforder = before surface tilt decenter order, 0 => dec then tilt, 1 =>
%     tilt then dec
% befdectilt = [dec-x dec-y tilt-x tilt-y tilt-z];
% aftstatus: 0 => explicit, 1=> pickup this surface, 2=> reverse this
%     surface, 3=> pickup surface # this-1, 4=> reverse surface # this-1,
%     etc.
% aftorder = after surface tilt decenter order, 0 => dec then tilt, 1 =>
%     tilt then dec
% aftdectilt = [dec-x dec-y tilt-x tilt-y tilt-z];
%
% strsurf = z_getsurftiltdec with no input arguments returns an empty
% struct

global MM P;

if nargin == 0,
    strsurf = struct('beforder',0,'befdectilt',zeros(1,5),'aftstatus',0,...
        'aftorder',0,'aftdectilt',zeros(1,5));
    return
end

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

DDEGETDATA = @(a)...
    (str2num(ddereq(zchan,['GetSurfaceData,' num2str(nsurf) ',' num2str(a)],[1 1])));

strsurf.beforder = DDEGETDATA(51);

for ii = 1:5,
    strsurf.befdectilt(ii) = DDEGETDATA(51+ii);
end
strsurf.befdectilt = strsurf.befdectilt .* [MM MM P P P];

strsurf.aftstatus = DDEGETDATA(60);
strsurf.aftorder = DDEGETDATA(61);

for ii = 1:5,
    strsurf.aftdectilt(ii) = DDEGETDATA(61+ii);
end
strsurf.aftdectilt = strsurf.aftdectilt .* [MM MM P P P];


return