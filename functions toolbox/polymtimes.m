function G = polymtimes(H1,H2)
% G = polymtimes(H1,H2)
% H1 and H2 are polynomial matrices (H(:,:,1) is first polynomial coeff., etc.)

N1 = size(H1); N2 = size(H2);
if any([length(N1) length(N2)] ~= [3 3]), error('H1 and H2 must be matrices of polynomials'); end
if N1(2) ~= N2(1), 
   error('number of columns of H1 must equal number of rows of H2');
end
Nr = N1(1);
Nc = N2(2);
Ni = N1(2);
Np = N1(3) + N2(3) - 1;

G = zeros(Nr,Nc,Np);

for n = 1:Nr,  % result row
   for m = 1:Nc,  % result column
      
      polytmp = zeros(Ni,Np);
      for ii = 1:Ni,
         polytmp(ii,:) = conv(H1(n,ii,:),H2(ii,m,:));
      end
      G(n,m,:) = sum(polytmp);
   end
end

return
      