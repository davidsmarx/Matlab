function [out1, out2] = z_getwave(zchan,nwave)
% [lambda,weight] = z_getwave(zchan,nwave)
% [primary, num_waves] = z_getwave(zchan,0)
%
% zchan from z_init
% nwave (optional) = integer wavelength number (default = 1)
% lambda = wavelength (m)
% weight = weight applied to this wavelength (default = 1.0)
% primary = nwave of primary wavelength
% num_waves = total number of wavelengths defined

global UM;

if nargin == 0, disp('usage: [lambda, weight] = z_getwave(zchan,nwave)'); return, end
if nargin == 1, nwave = 1; end

cmdstr = sprintf('GetWave,%d',nwave);
stmp = ddereq(zchan,cmdstr,[1 1]);
retval = sscanf(stmp,'%f,',[1 inf]);

if nwave == 0,
    primary = retval(1);
    Nwavelengths = retval(2);
    
    out1 = primary;
    out2 = Nwavelengths;
else,
    lambda = retval(1)*UM;
    weight = retval(2);
    
    out1 = lambda;
    out2 = weight;
end

return
