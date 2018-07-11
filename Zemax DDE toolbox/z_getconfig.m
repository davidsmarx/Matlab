function [cf, ncf] = z_getconfig(zchan)
% [cf, ncf] = z_getconfig(zchan)
% cf   = current configuration number
% ncf  = the number of configurations

try
cmdstr = ['GetConfig'];
a = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

cf = a(1);
ncf = a(2);

catch,
   disp('ERROR in z_getconfig!');
   keyboard;
end

return


