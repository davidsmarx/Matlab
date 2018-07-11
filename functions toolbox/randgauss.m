function y = randqauss(a)
% y = randgauss(a)
% returns an array of gaussian random numbers
% a is a size vector, see rand

x1 = rand(a);
x2 = rand(a);

y = sqrt(-2.0.*log(x1)).*cos(2.*pi.*x2);