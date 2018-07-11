function y = shift_array(x,n)
% y = shift_array(x)
% rotate the elements of x n places and wrap the last n elements to the front

nx = length(x);
y = x;
if n>0,
   y(1:n) = x(nx-n+1:nx);
   y(n+1:nx) = x(1:nx-n);
   
elseif n<0,
   m = -n;
   y(nx-m+1:nx) = x(1:m);
   y(1:nx-m) = x(m+1:nx);
end

