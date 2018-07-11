function y = besseljp(nu, r, k)
% y = besseljp(nu, r, k)
%
% besseljp is d besselj(nu,k*r) / dr

if nu == 0,
    y = -k.*besselj(1,k.*r);
    return
end

y = (nu./r).*besselj(nu,k.*r) - k.*besselj(nu+1, k.*r);

% resolve the r == 0 points
inan = abs(r) <= eps;

% for nu == 1, limit besseljp(1,0) => 1/2
% for nu > 1, limit besseljp(nu,0) => 0

if nu == 1,
    y(inan) = 0.5*k;
else
    y(inan) = 0.0;
end