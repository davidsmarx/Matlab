function y = mod2pi(x)
% y = mod2pi(x)
% y = x, -pi < x < pi
% y = x + 2k pi for some integer k, otherwise

y = x;
xp = (x > -pi)&(x < pi);

y(~xp) = mod( (x(~xp)+pi), 2*pi ) - pi;
