function opval = z_setNSCparm(zchan, nsurf, nobject, nparm, opval)
% opval = z_setNSCparm(zchan, nsurf, nobject, nparm, parmval)
%
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC
% nparm = parameter number

cmdstr = sprintf('SetNSCParameter,%d,%d,%d,%f',nsurf,nobject,nparm,opval);

opval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);


return
