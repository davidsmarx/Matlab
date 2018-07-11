function intensity = wavefront_intensity(A,dx,dy)
% intensity = wavefront_intensity(A,dx,dy)
% 
% integrates intensity of 2-d complex field A

[ny, nx] = size(A);

sumtmp = zeros(1,nx);

for ix = 1:nx,
    atmp = A(:,ix);
    sumtmp(ix) = atmp'*atmp;
end

intensity = sum(sumtmp)*dx*dy;
