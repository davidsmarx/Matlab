function coef = quadsolv(x1,x2,x3)

A = [x1.^2 x2.^2 x3.^2; x1 x2 x3; 1 1 1];
h = x3 - x1;
B = [h.*(x3.^2+x3.*x1+x1.^2)/3; h.*(x3+x1)/2; h];

coef = A\B;
