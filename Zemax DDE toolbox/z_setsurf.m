function errcond = z_setsurf(zchan, surfdata)
% errcond = z_setsurf(zchan, surfdata)
% surfdata is a struct with the following fields:
% nsurf, type, curv, thick, glass, radius, conic, parms, coat

global MM UM;

surfstr = ['SetSurface, ',...
   num2str(surfdata.nsurf) ',',...
   surfdata.type ', ',...
   num2str(surfdata.curv/(1/MM)) ', ',...
   num2str(surfdata.thick/MM) ',',...
   surfdata.glass ', ',...
   num2str(surfdata.radius/MM) ', ',...
   num2str(surfdata.conic) ', ',...
   sprintf('%e,',surfdata.parms),...
   surfdata.coat
];

retstr = ddereq(zchan,surfstr,[1 1]);

%disp(surfstr);
%disp(retstr);

% call GetUpdate to update solves, etc.
retval = ddereq(zchan,'GetUpdate',[1 1]);
if retval ~= '0', error('Lens Update failed!'); end

errcond = retval;
return