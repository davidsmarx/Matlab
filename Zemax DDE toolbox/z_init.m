function zchan = z_init
% zchan = z_init
%
% initialize DDE interface with ZEMAX
% and call GetRefresh

zchan = ddeinit('zemax','topic');
retstr = ddereq(zchan,'GetRefresh',[1 1]);

if sscanf(retstr,'%d') ~= 0,
   rc = ddeterm(zchan);
   error('Error: Update Failed!');
end

return