function [fm, xm, ym] = findpeak2(f01,f10,f11,f12,f21)
% [fm xm ym] = findpeak2(f01,f10,f11,f12,f21)
% [fm xm ym] = findpeak2(A)
%
% quadratic fit to interpolate peak value
% fnm are the values of some unknown function at evenly separated points
%          f01
%      f10 f11 f12
%          f21
% if input is a single matrix A, then findpeak2 fits the peak to the
% max2d(A) and the four pixels around it.
%
% fm is the value of the peak
% xm and ym is the location of the peak (w.r.t. the location of f11)
% f11 must be max(fnm) to find peak or
% f11 must be min(fnm) to find minimum

% AI is the inverse of the matrix A: Anm = [xn^2 xn yn^2 yn 1]
% and fnm = A * [a;b;c;d;e], based on the parabaloid equation:
% f = a x^2 + b x + c y^2 + d y + e

if nargin == 1,
    [Nr, Nc] = size(f01);
    [amax, nr, nc] = max2d(f01);
    if any([nc==1, nc==Nc, nr==1, nr==Nr]),
        warning('peak occurs at edge of matrix');
        [fm, xm, ym] = deal(NaN);
        return
    end
    vp = [f01(nr-1,nc);f01(nr,nc-1);f01(nr,nc);f01(nr,nc+1);f01(nr+1,nc)];
elseif nargin == 5,
    nr = 0; nc = 0;
    vp = [f01;f10;f11;f12;f21];
else,
    error('invalid input');
end

AI = [0 0.5 -1 0.5 0;0 -0.5 0 0.5 0;0.5 0 -1.0 0 0.5;-0.5 0 0 0 0.5;0 0 1 0 0];
a  = AI * vp;

xmp = -0.5*a(2)/a(1);
ymp = -0.5*a(4)/a(3);
fm = [xmp.^2 xmp ymp.^2 ymp 1] * a;

xm = nc + xmp;
ym = nr + ymp;

return
