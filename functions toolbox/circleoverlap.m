function overlaparea = circleoverlap(r0,R1,R2)
% overlaparea = circleoverlap(r0,R1,R2)
%
% calculate the overlap area of two circles, see p. 43 of notebook
% r0 = distance between centers of the two circles
% R1 = radius of larger circle
% R2 = radius of smaller circle

% 2005-11-14

if R2 > R1, rtmp = R2; R2 = R1; R1 = rtmp; end

t2 = @(r0)(acos( (R1.^2 - R2.^2 - r0.^2)./(2.*r0.*R2) ) );
t1 = @(r0)((1/2).*acos( (R2.^2./R1.^2) .* (cos(2*t2(r0)) - 1) + 1) );

A1 = @(r0)(R1.^2 .* (t1(r0) - sin(2*t1(r0))/2) );
A2 = @(r0)(R2.^2 .* (pi - t2(r0) + sin(2*t2(r0))/2) );

overlaparea = zeros(size(r0));
overlaparea(r0 <= R1 - R2) = pi*R2.^2;
ib = r0 > R1 - R2 & r0 <= R1 + R2;
overlaparea(ib) = A1(r0(ib)) + A2(r0(ib));

end
