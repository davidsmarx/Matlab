function mtf = MTFtheory(ff,fcutoff)
% mtf = MTFtheory(ff,fcutoff)
%
% equation 10, Vol. 2, Ch. 32, Handbook of Optics

switch nargin,
    case 0,
        ff = linspace(0,1)';
        fcutoff = 1;
    case 1,
        fcutoff = 1;
    case 2,
    otherwise,
end

k = abs(ff)./fcutoff;

n0 = k<=1;
k0 = k(n0);
n1 = k>1;

mtf = zeros(size(k));

mtf(n0) = (2./pi).*( acos(k0) - k0.*sqrt(1-k0.^2) );
mtf(n1) = 0;

return