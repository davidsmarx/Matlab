function m = maxabs(x)
% m = maxabs(x)
%
% m = the element of x with the largest absolute value
% 
% [xmax, imax] = max(abs(x(:)));
% m = sign(x(imax)) * xmax;


[xmax, imax] = max(abs(x(:)));
m = sign(x(imax)) * xmax;
