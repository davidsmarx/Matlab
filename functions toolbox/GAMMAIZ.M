% GAMMAIZ incomplete gamma function that also allows for complex "a" 
%        arguments, unlike MATLAB's GAMMAINC function. For real arguments, 
%        GAMMAIZ gives the same results as GAMMAINC (within numerical
%        error).
%
%   usage: f = gammainc(x,a)      x, a can be matrices of same size or scalar
%   
%   algorithm: slight modification of MATLAB's GAMMAINC: call GAMMALNZ
%     instead of GAMMALN

% NOTICE swapped arguments in GAMMAINC(x,a), whereas custom notation is 
% "P(a,x)".

function b = gammainc(x,a)

if nargin~=2, error('Requires two input arguments.'); end
if prod(size(x))==1, x = x(ones(size(a))); end
if prod(size(a))==1, a = a(ones(size(x))); end

b = x;
lngama = gammalnz(a);

k = find(x == 0);
if ~isempty(k)
   b(k) = 0*k;
end

k = find((x ~= 0) & (x < a+1));
if ~isempty(k)
   if prod(size(a)) == 1
      % Series expansion for x < a+1
      ap = a;
      sum = ones(size(k))/a;
      del = sum;
      while norm(del,'inf') >= 10*eps*norm(sum,'inf')
         ap = ap + 1;
         del = x(k) .* del/ap;
         sum = sum + del;
      end
      b(k) = sum .* exp(-x(k) + a*log(x(k)) - lngama);
   else
      % Series expansion for x < a+1
      ap = a;
      sum = ones(size(k)) ./a(k);
      del = sum;
      while norm(del,'inf') >= 10*eps*norm(sum,'inf')
         ap = ap + 1;
         del = x(k) .* del./ap(k);
         sum = sum + del;
      end
      b(k) = sum .* exp(-x(k) + a(k) .*log(x(k)) - gammalnz(a(k)));
   end
end

k = find(x >= a + 1);
if ~isempty(k)
   % Continued fraction for x >= a+1
   a0 = ones(size(k));
   a1 = x(k);
   b0 = zeros(size(k));
   b1 = a0;
   fac = 1;
   n = 1;
   g = b1;
   gold = b0;
   if prod(size(a)) == 1
      while norm(g-gold,'inf') >= 10*eps*norm(g,'inf');
         gold = g;
         ana = n - a;
         a0 = (a1 + a0*ana) .* fac;
         b0 = (b1 + b0*ana) .* fac;
         anf = n*fac;
         a1 = x(k) .* a0 + anf .* a1;
         b1 = x(k) .* b0 + anf .* b1;
         fac = 1 ./ a1;
         g = b1 .* fac;
         n = n+1;
         b(k) = 1 - exp(-x(k) + a*log(x(k)) - lngama) .* g;
      end
   else
      while norm(g-gold,'inf') >= 10*eps*norm(g,'inf');
        gold = g;
        ana = n - a(k);
        a0 = (a1 + a0 .*ana) .* fac;
        b0 = (b1 + b0 .*ana) .* fac;
        anf = n*fac;
        a1 = x(k) .* a0 + anf .* a1;
        b1 = x(k) .* b0 + anf .* b1;
        fac = 1 ./ a1;
        g = b1 .* fac;
        n = n+1;
        b(k) = 1 - exp(-x(k) + a(k) .* log(x(k)) - gammalnz(a(k))) .* g;
      end
   end
end

