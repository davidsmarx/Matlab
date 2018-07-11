% GAMMALNZ log gamma function that also allows for complex arguments, unlike
%        MATLAB's LOGGAMMA function. For real arguments, GAMMALNZ gives the
%        same results as LOGGAMMA (within numerical error).
%
%   algorithm: numerical recipes 6.1.5, "gammln"
%
%   usage: [f] = gammaz(z)    z can be any size

function [f] = gammalnz(z);

[m n] = size(z);
f = zeros(m,n);
cof = [76.18009172947146
      -86.50532032941677
       24.01409824083091
       -1.231739572450155
        0.1208650973866179e-2
       -0.5395239384953e-5];		% 6 coefficients in series expansion

for icol = 1:n,				% do one column at a time
  zz = z(:,icol);
  zp = ones(6,1)*zz' + [1:length(cof)]'*ones(1,length(zz)); % vectorize
  ser = (cof'*(1./zp)+1.000000000190015)';
  tmp = zz+5.5 - (zz+.5).*log(zz+5.5);
  f(:,icol) = -tmp + log(2.5066282746310005*ser./zz);
end

   


