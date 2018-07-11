function y = rect(x)
% y = rect(x)
% y = 1, x is in (-1/2,1/2]
% y = 0, otherwise

y = (x>-0.5) & (x<=0.5);