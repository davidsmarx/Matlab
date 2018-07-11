function otf = OTF_incoherent_circle(r)
% otf = OTF_incoherent_circle(r)
%
% r is normalized aperture radius
%
% equation 6-31 on p. 120 of Goodman

z = r<=2;

otf = zeros(size(r));

r2 = 0.5*r(z);

otf(z) = (2./pi).*( acos(r2) - r2.*sqrt(1 - r2.^2));
