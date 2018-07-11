function opval = z_getNSCparm(zchan, nsurf, nobject, nparm)
% opval = z_getNSCparm(zchan, nsurf, nobject, nparm
%
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC
% nparm = parameter number


cmdstr = sprintf('GetNSCParameter,%d,%d,%d',nsurf,nobject,nparm);
opvalstr = ddereq(zchan,cmdstr,[1 1]);

opval = sscanf(opvalstr,'%f,',[1 inf]);

return
