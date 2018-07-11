function [U, Y] = thinlensraytrace(d,p,u0,y0)
% [U, Y] = thinlensraytrace(d,p,u0,y0)
%   (used to be 'paraxraytrace_lens()')
%
% d = list of lens separations
% p = list of lens powers (p = 1/f)
% length(d) must equal length(p). d(n) = distance from p(n) to p(n+1)
% u0,y0 = slope and height at entrance to first lens
% at each lens, apply lens then distance to next lens
% U, Y = vectors of slope, height at entrance to each lens

if length(d) ~= length(p), error('d and p must be same length'); end

A = sym([u0; y0]);

for i = 1:length(d),
   A(:,end+1) = [1 -p(i);d(i) 1-d(i)*p(i)]*A(:,end);
end

U = A(1,:);
Y = A(2,:);

return