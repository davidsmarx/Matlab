function [R, T] = TFfilter_response(n,d,theta,lam,tetm)
% [R, T] = TFfilter_response(n,d,theta,lam,tetm)
% R, T are complex reflected and transmitted field coefficients
% size(R) = size(T) = [length(lam), length(theta)]
% n(1) = index of substrate
% n(end) = index of incident medium
% length(d) = length(n)-2
% theta [rad]
% tetm: 0 => TE (default), 1 => TM, (or 'S' = TE, 'P' = TM)

if nargin == 0, error('usage: [R, T] = TFfilter_response(n,d,theta,lam,tetm)'); end

if ~exist('tetm'), 
   tetm = 0; 
else,
   if ischar(tetm),
      if strcmp(tetm,'S'), tetm = 0;
      elseif strcmp(tetm,'P'), tetm = 1;
      else, error(['unknown polarization: ' tetm]);
      end
   end
end

theta = theta(:);
lam = lam(:);

R = zeros(length(lam),length(theta));
T = R;

n = flipud(n(:)); % thin_film_filter_2 uses opposite convention, n(1) = incident
d = flipud(d(:));

for it = 1:length(theta),
   
   for il = 1:length(lam),
      [r t rr tt] = thin_film_filter_2(n,d,theta(it),lam(il),tetm);
      R(il,it) = rr(1);
      T(il,it) = tt(end);
   end
   
end

return