function [U, Y, I] = paraxraytrace_surf(d,c,n,A1)
% [U, Y, I] = paraxraytrace_surf(d,c,n,A1)
% d = list of surface separations
% c = list of surface curvatures (c = 1/R)
% n = list of medium indexes
% length(d) = # surfaces-1
% length(c) = # surfacas
% length(n) = # surfaces+1
% A1 = [u0;y0] = slope and height at entrance to first lens
%
% at each surface, apply surface then distance to next surface
% U, Y = vectors of slope, height at exit of each surface
% I = angle of incidence into each surface


%p = c.*diff(n)./n(2:end); % calculate power of each surface
p = c.*(n(2:end)-n(1:end-1))./n(2:end);
% set first thickness = 0 (starting in front of first surface)
d = [0; d(:)];

% % imitialize A (reverse convention)
if isa([d(:); c(:); n(:); A1],'sym'),
   A = sym([A1(2); A1(1)]);
else,
   A = [A1(2); A1(1)];
end

for i = 1:length(d),
   A(:,end+1) = [1 d(i);-p(i) n(i)./n(i+1)-d(i)*p(i)]*A(:,i);
   I(i) = c(i)*A(1,i+1) + A(2,i); % y after surface, u before surface
end

Y = A(1,2:end);
U = A(2,2:end);

if isa([Y; U],'sym'),
   Y = simplify(Y);
   U = simplify(U);
end

return