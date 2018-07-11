function errcond = z_setsurfparm(zchan, surfdata, nparm)
% errcond = z_setsurfparm(zchan, surfdata, nparm)
% surfdata is a struct with the following fields:
% nsurf, type, curv, thick, glass, radius, conic, parms, coat
% new parameter data is in surfdata.parms
% nparm is the parameter number to be set

surfstr = ['SetSurfaceParameter, ',...
   num2str(surfdata.nsurf) ',',...
   num2str(nparm) ', ',...
   num2str(surfdata.parms(nparm))
];

retstr = ddereq(zchan,surfstr,[1 1]);

%disp(surfstr);
%disp(retstr);

errcond = retstr;
return