function opval = z_setNSCproperty(zchan, nsurf, nobject, property_code, face, property_string)
% opval = z_setNSCproperty(zchan, nsurf, nobject, property_code, face, property_string)
%
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC
% code = specifies the property, see p.790, "SETNSCPROPERTY" in the ZPL chapter
% face = # from the drop down list of faces
%
% SetNSCProperty, surface, object, code, face, value

cmdstr = sprintf('SetNSCProperty,%d,%d,%d,%d,%s',...
    nsurf,nobject,property_code,face,property_string);

opval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);


return
