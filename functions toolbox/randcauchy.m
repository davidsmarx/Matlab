function y = randcauchy(a)
% y = randcauchy(a)
% returns an array of cauchy distributed random numbers
% f(y) = pi/(1 + y^2)
% a is a size vector [ncolumns nrows]

phi = (pi/2)*(rand(a) - 0.5);
y = tan(phi);

return