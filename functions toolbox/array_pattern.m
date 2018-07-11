function f = array_pattern(x,N)
% f = array_pattern(x,N)
% returns the array pattern function sin(pi N x)/sin(pi x)

if length(N)~=1,
   error('N must be a scalar')
end
if N<1,
   error('N must be > 1')
end
if mod(N,1)>1.0e-6,
   error('N must be an integer')
end


f = zeros(size(x));

i_even = abs(mod(x,2))<1.0e-9;
i_odd  = abs(mod(x+1,2))<1.0e-9;
f(i_even) = N;
if mod(N,2)<1.0e-6, %N is even
   f(i_odd) =-N;
else
   f(i_odd) = N;
end

i_else = ~(i_even | i_odd);
f(i_else) = sin(pi*N*x(i_else)) ./ sin(pi*x(i_else));

return
