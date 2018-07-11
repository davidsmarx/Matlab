function status = z_getupdate(zchan)
% status = z_getupdate(zchan)
%
% status: 0 => success, -1 => failure

if nargin == 0, disp('usage: status = z_getupdate(zchan)'); return, end

% call update
retstr = ddereq(zchan,'GetUpdate',[1 1]);

status = str2num(retstr);

return
