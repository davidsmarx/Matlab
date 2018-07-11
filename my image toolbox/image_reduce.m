function B = image_reduce(A,P,Q)
% B = image_reduce(A,P,Q)
%
% resample input image A at rate P/Q

[Ny, Nx, N3] = size(A);

Px = P; Qx = Q;
Py = P; Qy = Q;

% use interp2
xi = [1:Qx/Px:Nx];
yi = [1:Qy/Py:Ny]';
for ii = 1:N3,
   dtmp(:,:,ii) = interp2(double(A(:,:,ii)),xi,yi);
end

switch class(A),
case 'uint8',
   B = uint8(round(dtmp));
case 'double',
   B = dtmp;
otherwise,
   error(['not compatible with ' class(A)]);
end

return