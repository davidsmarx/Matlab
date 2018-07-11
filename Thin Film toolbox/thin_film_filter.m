function [R, T, rr, tt] = thin_film_filter(n0,n,d,theta,lam,tetm)
% [R, T] = thin_film_filter(n0,n,d,theta,lam)
% n0 = index of incident medium
% n = index of refraction for each layer. 
%     n(N) = index of transmission medium
% d = thickness of each layer, not counting incident medium or transmission
%     medium. length(d) = length(n)-1
%     So, n(1) d(1) are the index and thickness of the first thin film,
%     n(N-1) d(N-1) are the index and thickness of the last film,
%     and n(N) is the index of the transmission medium
% theta = angle of incidence [rad]
% lam = wavelength. units of lam must be same as d
% R = vector of complex field reflection coefficients at each layer
% T =          "              transmission        "
% tetm: 0 => TE (default), 1 => TM
% R = vector of normalized reflected intensity coefficient
% T =        "             transmitted     "
% rr = vector of complex field reflection coefficient
% tt =      "                  transmission "

if nargin < 5, tetm = 0; end

N = length(n);
if length(d) ~= N-1, error('n and d mismatch'); end
n = n(:); d = d(:);

kx = 2*pi*n0*sin(theta)/lam;
kz0 = sqrt( (2*pi*n0/lam).^2 - kx.^2);
kz = sqrt( (2*pi*n/lam).^2 - kx.^2 ); 
if tetm == 1,
   kzz0 = kz0./(n0.^2);
   kzz = kz./(n.^2);
else,
   kzz0 = kz0;
   kzz = kz;
end

eep = exp(j*kz(1:N-1).*d);
eem = exp(-j*kz(1:N-1).*d);

nnz = 8*N-4;  % number of non-zero elements in the sparse matrix
i = zeros(nnz,1);  % row index
j = zeros(nnz,1);  % column index
s = zeros(nnz,1);  % sparse matrix elements
i(1:6) = [ones(1,3), 2*ones(1,3)].';
j(1:6) = [1:3 1:3].';
s(1:6) = [1 -1 -1 kzz0 kzz(1) -kzz(1)].';
for ll = 2:N-1,
   ii = 6+8*(ll-2);
   i(ii+1:ii+8) = [(2*ll-1)*ones(1,4), (2*ll)*ones(1,4)].';
   j(ii+1:ii+8) = [2*ll-3 + [1:4], 2*ll-3 + [1:4]].';
   s(ii+1:ii+8) = [ eem(ll-1) eep(ll-1) -1 -1 -kzz(ll-1)*eem(ll-1) kzz(ll-1)*eep(ll-1) kzz(ll) -kzz(ll) ].';
end
i(nnz-5:nnz) = [(2*N-1)*ones(1,3), (2*N)*ones(1,3)].';
j(nnz-5:nnz) = [2*N-3 + [1:3], 2*N-3 + [1:3]].';
s(nnz-5:nnz) = [eem(N-1) eep(N-1) -1 -kzz(N-1)*eem(N-1) kzz(N-1)*eep(N-1) kzz(N)].';

S = sparse(i,j,s);
b = [-1 kzz0 zeros(1,2*N-2)]';

x = S \ b;
rr = x(1:2:2*N);
tt = x(2:2:2*N);
R = abs(rr).^2;
if tetm == 1,
   R = [abs(rr(1)).^2; real( (kz(1:N-1)./(n(1:N-1).^2)) ./ (kz0./(n0.^2)) ) .* (abs(rr(2:end)).^2)];
   T = real( (kz./(n.^2)) ./ (kz0./(n0.^2)) ) .* (abs(tt).^2);
else,
   R = [abs(rr(1)).^2; real( (kz(1:N-1)./kz0) ) .* (abs(rr(2:end)).^2)];
   T = real( (kz./kz0) ) .* (abs(tt).^2);
end

% check that:
% Rtm(1) + Ttm - [Rtm(2:end); 0] = ones(N,1)

return

if 0,  % the same routine without using sparse matrix
   
A = zeros(2*N);
A(1,1:3) = [1 -1 -1];
A(2,1:3) = [kz0 kz(1) -kz(1)];
for i = 2:N-1,
   A(2*i-1,2*i-2:2*i+1) = [ eem(i-1) eep(i-1) -1 -1 ];
   A(2*i,2*i-2:2*i+1)   = [-kz(i-1)*eem(i-1) kz(i-1)*eep(i-1) kz(i) -kz(i) ];
end
A(2*N-1,2*N-2:2*N) = [ eem(N-1) eep(N-1) -1 ];
A(2*N,2*N-2:2*N)   = [-kz(N-1)*eem(N-1) kz(N-1)*eep(N-1) kz(N) ];

b = [-1 kz0 zeros(1,2*N-2)]';

rr = A \ b; % solves matrix equation A * rr = b

R = rr(1:2:2*N);
T = rr(2:2:2*N);

end