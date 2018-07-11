function [fm, xm] = findpeak(f1,f2,f3)
% [fm xm] = findpeak(f1,f2,f3), f1, f2, f3 are vectors of the same size
% [fm xm] = findpeak(ff), where size(ff) = [3 N]
% quadratic fit to interpolate peak value
% f1, f2, f3 are the values of some unknown function at evenly separated points
% fm is the value of the peak
% xm is the location of the peak (w.r.t. the location of f2)
% f2 must be > than both f1 and f3 to find peak or
% f2 must be < than both f1 and f3 to find minimum

if nargin == 0, error('usage: [fm, xm] = findpeak(f1,f2,f3)'); end

if nargin == 1,
   f3 = f1(3,:); f2 = f1(2,:); f1 = f1(1,:);
end
   
a = 0.5*(f1+f3) - f2;
b = 0.5*(f3-f1);
c = f2;

xm = -0.5*b ./ a;
fm = a .* (xm.^2) + b .* xm + c;

return
