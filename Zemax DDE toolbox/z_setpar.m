function er = z_setpar(zchan,surf,parmlist,vallist)
% status = z_setpar(zchan,surf,parmlist,vallist)
% surf = surface struct: 
%          nsurf, type, curv, thick, glass, radius, conic, parms, coat

surf.parms(parmlist) = vallist;
er = z_setsurf(zchan, surf);

return

