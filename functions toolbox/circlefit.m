function [r, x0, y0] = circlefit(x,y)
% [r, x0, y0] = circlefit(x,y)
%
% least squares solution to find the center (x0,y0) and radius r that
% minimizes the difference between the distance from given points x,y to
% (x0,y0) and the radius r

uu = [x(:) y(:) ones(size(x(:)))] \ (0.5*[x(:).^2 + y(:).^2]);
x0 = uu(1);
y0 = uu(2);
r  = sqrt(2*uu(3) + x0.^2 + y0.^2);
