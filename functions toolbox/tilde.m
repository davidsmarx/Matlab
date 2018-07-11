function ht = tilde(h)
% ht = tilde(h)
% h is a polynomial or matrix of polynomials (FIR only)
% tilde(h) = fliplr(conj(h))

N = size(h);

switch length(N)
case 2,
   if length(h)==length(h(:)), % must be a single polynomial
      ht = fliplr(h(:)');   % output is a row vector polynomial
   else, % h is a vector of polynomials
      for n = 1:N(1),
         ht(n,:) = fliplr(conj(h(n,:)));
      end
   end
case 3,
   for n = 1:N(3),
      ht(:,:,n) = h(:,:,n)';
   end
   ht = flipdim(ht,3);
otherwise,
   error(['tilde for ' num2str(N) ' is undefined']);
end

return

