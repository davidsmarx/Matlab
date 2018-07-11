function otf = OTF_partiallycoherent_circle(R1,R2,r)
% otf = OTF_partiallycoherent_circle(R1,R2,r)
%
% R1 = radius of larger aperture (typically objective lens pupil)
% R2 = radius of smaller aperture (typically aperture stop)
%
% p. 43 of notebook

otf = zeros(size(r));

% region where smaller circle is entirely within larger circle
n1 = r <= R1 - R2;
otf(n1) = 1;

% cross-over region
n1 = r > R1 - R2 & r < R1 + R2;
cost2 = (R1.^2 - R2.^2 - r(n1).^2)./(2*r(n1)*R2);
if any(abs(cost2)>1), error('t2'); end
t2 = acos(cost2);

cos2t1 = (R2./R1).^2 * (cos(2*t2)-1) + 1;
if any(abs(cos2t1)>1), error('t1'); end
t1 = acos(cos2t1)/2;

A1 = R1.^2*(t1 - 0.5*sin(2*t1));
A2 = R2.^2*(pi - t2 + 0.5*sin(2*t2));
otf(n1) = (A1 + A2)./(pi*R2.^2);

% outer region is everything else

return