function [val,bp,wf]=gaussint(fun,a,b,n)
% [val,bp,wf]=gaussint(fun,a,b,n) integrates
% a function from a to b using an n-point
% Gauss rule which is exact for a polynomial
% of degree 2*n-1. Concepts on page 93 of
% 'Methods of Numerical Integration' by
% Philip Davis and Philip Rabinowitz yield
% the base points and weight factors.
%
% 'fun' -  string variable for the name of the function being integrated
% a,b   -  integration limits
% n     -  order of formula (default is 20)
% val   -  value of the integral
% bp,wf -  Gauss base points and weight
%          factors
if nargin == 0,
   fprintf('usage: [val,bp,wf] = gaussint(fun,a,b,n)\n');
   return
end
 
if nargin < 4, n=20; end
u=(1:n-1)./sqrt((2*(1:n-1)).^2-1);
[vc,bp]=eig(diag(u,-1)+diag(u,1));
[bp,k]=sort(diag(bp));
wf=2*vc(1,k)'.^2;
x=(a+b)/2+(b-a)/2*bp;
f=feval(fun,x)*(b-a)/2; val=wf(:)'*f(:);