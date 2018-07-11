function retval = z_InsertObject(zchan,nsurf,nobject)
% retval = z_InsertObject(zchan,nsurf,nobject)
%
% nsurf: surface
% nobject: object number
%
% retval = 0 if success, NaN if failed

cmdstr = ['InsertObject,' num2str(nsurf) ',' num2str(nobject)];

retstr = ddereq(zchan,cmdstr,[1 1]);
retval = str2double(retstr);

return


