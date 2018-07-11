function phi = interpphase(x1,phi1,xi)
% phi = interpphase(x1,phi2,xi)
% linear interpolation using mod2pi so that the
% interpolated value is always on the shortest branch.
% phi1 are given values at points x1.
% xi are the points at which to interpolate phi1
% each column of phi1 is interpolated separately

if nargin ~= 3,
   error('usage: phi = interpphase(x1,phi1,xi)');
end

if ~isreal(phi1) | ~isreal(x1) | ~isreal(xi),
   error('input phases and abscissae must be real');
end

[Nx1 Nc] = size(phi1);
if length(x1)~=Nx1, error('length(xi) must equal # rows of phi1'); end
Nxi = length(xi);

phi = zeros(Nxi,Nc);
for i=1:Nxi,
   
   [dmin, imin] = min(abs(xi(i)-x1));
   
   % nearest neighbor:
   phi(i,:) = phi1(imin,:);
   %

   if xi(i)-x1(imin) < 0,
      ia = imin-1;
      ib = imin;
   else
      ia = imin;
      ib = imin+1;
   end
   if (ia<1) | (ib>Nx1), 
      %error('interpolation index out of bound');
      %fprintf('index at boundary: ia=%d ib=%d\n',ia,ib);
      if ia<1, phi(i,:) = phi1(ib,:);
      else phi(i,:) = phi1(ia,:);
      end
   else
      
      if ia==ib, error('ia = ib'); end
      xa = x1(ia);
      xb = x1(ib);
      if (xi(i)-xa)<0, 
         warning('x1 is not monotonic'); 
         xa = xi(i);
      end
      if abs(xa-xb)<1.0e-12, 
         fprintf('ia = %d ib = %d\n',ia,ib);
         fprintf('x1 = %f x1 = %f\n',x1(ia),x1(ib));
         error('delta too small'); 
      end
      xx = (xi(i)-xa)/(xb-xa);
      phi(i,:) = mod2pi(phi1(ia,:) + xx .* mod2pi(phi1(ib,:) - phi1(ia,:)));
   end
   
end


%phi = mod2pi( phi1 + x .* mod2pi(phi2 - phi1) );
