function [R, T] = fitler_response(n,d,theta,lam,tetm)
% [R, T] = fitler_response(n,d,theta,lam,tetm)
% R, T are complex reflected and transmitted field coefficients
% size(R) = size(T) = [length(lam), length(theta)]

if ~exist('tetm'), tetm = 0; end

theta = theta(:);
lam = lam(:);

R = zeros(length(lam),length(theta));
T = R;

for it = 1:length(theta),
   
   for il = 1:length(lam),
      [r t rr tt] = thin_film_filter_2(n,d,theta(it),lam(il),tetm);
      R(il,it) = rr(1);
      T(il,it) = tt(end);
   end
   
end

return