function F = dct_(f,isign)
% F = dct_(f,isign)
% assumes # of rows(f) = 2^N
% isign = 1 for forward transform
% isign = -1 for inverse transform

[nr nc] = size(f);
F = zeros([nr nc]);

[j k] = meshgrid(0:nr-1,0:nr-1);
if isign == 1
	for i=1:nc
   	F(:,i) = cos((0.5*pi/nr)*((2*j+1).*k))*f(:,i);
   end
elseif isign == -1
   for i=1:nc
      F(:,i) = (2/nr)*(cos((0.5*pi/nr)*(j.*(2*k+1)))*f(:,i) - 0.5*f(1,i));
   end
else
   error('isign must be 1 or -1')
end

