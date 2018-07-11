function status = z_setsurftiltdec(zchan,nsurf,strsurf)
% status = z_setsurftiltdec(zchan,nsurf,structsurf)
% 
% nsurf is the surface #
% structsurf is a surface struct as returned by z_getsurftiltdec
%
% structsurf is a struct with the following fields:
% beforder = before surface tilt decenter order, 0 => dec then tilt, 1 =>
%     tilt then dec
% befdectilt = [dec-x dec-y tilt-x tilt-y tilt-z];
% aftstatus: 0 => explicit, 1=> pickup this surface, 2=> reverse this
%     surface, 3=> pickup surface # this-1, 4=> reverse surface # this-1,
%     etc.
% aftorder = after surface tilt decenter order, 0 => dec then tilt, 1 =>
%     tilt then dec
% aftdectilt = [dec-x dec-y tilt-x tilt-y tilt-z];

% set units
global MM P;
strsurf.befdectilt = strsurf.befdectilt./[MM MM P P P];
strsurf.aftdectilt = strsurf.aftdectilt./[MM MM P P P];

DDESETDATA = @(a,b)...
    (ddereq(zchan,['SetSurfaceData,' num2str(nsurf) ',' num2str(a) ',' num2str(b)],[1 1]));

status = DDESETDATA(51,strsurf.beforder);

for ii = 1:5,
    status = DDESETDATA(51+ii,strsurf.befdectilt(ii));
end

status = DDESETDATA(60,strsurf.aftstatus);
status = DDESETDATA(61,strsurf.aftorder);

for ii = 1:5,
    status = DDESETDATA(61+ii,strsurf.aftdectilt(ii));
end

return