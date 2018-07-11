function [gd, cd] = group_delay(R,df,f0)
% [gd, cd] = group_delay(R,df,f0)
%
% R  = complex sampled signal
% df = sample spacing (default = 1)
% f0 = center frequency (optional, used to calculate cd)
% gd = gradient(phase(R))/2pi [1/units(df)] (if df = [THz], then gd = [ps])
% cd = 1e-3*gradient(gd), 
%      if df = [THz] and f0 = [THz], then cd = [ps/nm]
%      if df = [THz] and no f0 is input, then cd = [ps/THz]

try df; catch df = 1; end

constants;

% unwrap phase from complex signal
r_phi = unwrap(angle(R));
r_phi = r_phi - mean(r_phi);

% gradient for GD
gd = (0.5/pi)*gradient(r_phi,df);

% gradient again for CD
try,
   N = length(R);
   ff = f0 + df*[-floor(N/2):ceil(N/2-1)]';
   
   cd = gradient(gd,df).*(-(ff.^2)./(1000*C));
catch,
   cd = gradient(gd,df);
end

return
