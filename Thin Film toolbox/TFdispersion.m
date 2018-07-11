function [gd, chrdisp, frqdisp] = TFdispersion(tt,x,xtype)
% [gd, chrdisp, frqdisp] = TFdispersion(tt,x,xtype)
%
% tt = complex filter response
% xtype = 'freq' or 'lam'
%
% gd is [ps]
% chrdisp is [ps/nm] if freq is [THz] or lam is [um]
% frqdisp is [ps/GHz]
if nargin == 0, error('usage: [gd, chrdisp, frqdisp] = TFdispersion(tt,x,xtype)'); end
global C;

y1 = C./x(end);
y2 = C./x(1);
y  = linspace(y1,y2,length(x));

if strcmp(xtype,'freq'),
	freq = x;
	lam  = y;	
	ttfr = tt;
elseif strcmp(xtype,'lam'),
	lam  = x;
	freq = y;
	ttfr = interp1(C./lam,tt,freq,'spline');
else,
	error(['unknown xtype: ' xtype]);
end
df = abs(freq(2)-freq(1));

% group delay
t_phi = unwrap(angle(ttfr));
t_phi = t_phi - mean(t_phi);
gd    = (0.5/pi)*gradient(t_phi,freq(2)-freq(1)); % [ps]

% frequency dispersion
frqdisp = 0.001*gradient(gd,df); % [ps/GHz]

% chromatic dispersion
c_gd = interp1(freq,gd,C./lam,'spline');
chrdisp = 0.001*gradient(c_gd,abs(lam(2)-lam(1))); % [ps/nm]

if any(isnan(gd)), disp('NaN in gd'); keyboard; end
if any(isnan(chrdisp)), disp('NaN in chrdisp'); keyboard; end
if any(isnan(frqdisp)), disp('NaN in frqdisp'); keyboard; end

return