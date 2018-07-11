function s = chebyshevpuls(t,T,a)
% s = chebyshevpuls(t,T,a)
% return a pulse formed from a chebyshev window with width T
% and sidelobe level -a dB. The default value for a is 30 dB
% t = time is a real 1-d vector
% if t is a matrix, return the chebyshevpuls of each column

if nargin < 2,
   error('usage: s = chebyshevpuls(t,T,a)')
elseif nargin == 2,
   a = 30;
end

[nr nc] = size(t);

if nr*nc == length(t),  % t must be a 1-d vector
   s = chebyshevpuls_vector(t,T,a);
   
else
   s = zeros([nr nc]);
   for i = 1:nc, % do it for each column
      s(:,i) = chebyshevpuls_vector(t(:,i),T,a);
   end
   
end % if t is a vector

      
function y = chebyshevpuls_vector(t,T,a)

y = 0 * t;

tp = (t>=-T/2)&(t<T/2);
y(tp) = chebwin(length(t(tp)),a);
