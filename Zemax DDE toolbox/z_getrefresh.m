function status = z_getrefresh(zchan)
% status = z_getrefresh(zchan)
%
% status: 0 => success, -1 => failure
%
% GetRefresh causes Zemax to copy the lens data from the LDE into the
% stored copy of the server.

if nargin == 0, disp('usage: status = z_getrefresh(zchan)'); return, end

% call refresh
retstr = ddereq(zchan,'GetRefresh',[1 1]);

status = str2num(retstr);

return
