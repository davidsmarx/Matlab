function [status, numberconfigs] = z_setconfig(zchan,nconfig)
% [status, numberconfigs] = z_setconfig(zchan,nconfig)
if nargin == 0, disp('usage: status = z_setconfig(zchan,nconfig)'); return, end

cmdstr = sprintf('SetConfig,%d',nconfig);
retval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

% retval(1) = nconfig
% retval(2) = number of configurations
% retval(3) = error code (0 = no error, -1 = cannot trace a ray)

if retval(1) ~= nconfig, error('error setting config'); end
if retval(3) ~= 0, error('SetConfig returned error'); end

numberconfigs = retval(2);
status = retval(3);

return
