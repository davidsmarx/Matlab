
function [dm1 dm2]	= get_dm_volt(fn);
%dd   = fitsinfo(fn);
%ndm  = length(dd.Image); 
dm1      = fitsread(fn,'image',1); 
dm1      = dm1(:,:,1);
dm2      = fitsread(fn,'image',2); 
dm2      = dm2(:,:,1);
return
