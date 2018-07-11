function y = markov(n_r,n_c,sig)
% y = markov(n_r,n_c,sig)
% y is a matrix where each column is a separate markov sequence

if nargin==0,
   disp('usage: y = markov(n_r,n_c,sig)')
   return
end

y = zeros(n_r, n_c);
r = [1 n_c];

y(1,:) = randgauss(r);
for i=2:n_r,
   y(i,:) = y(i-1,:) + randgauss(r);
end

