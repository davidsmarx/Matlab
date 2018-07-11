function opvalstr = z_getNSCproperty(zchan, nsurf, nobject, property_code, face)
% opval = z_getNSCparm(zchan, nsurf, nobject, nparm
%
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC
% code = specifies the property, see p.790, "SETNSCPROPERTY" in the ZPL chapter
% face = # from the drop down list of faces
%
% GetNSCProperty, surface, object, code, face
% 
% return value:
%     a string
%     ascii(10) for invalid face
%     '0.0000E00' for invalid property_code



cmdstr = sprintf('GetNSCProperty,%d,%d,%d,%d',nsurf,nobject,property_code,face);
opvalstr = ddereq(zchan,cmdstr,[1 1]);

return
