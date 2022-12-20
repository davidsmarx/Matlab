function Im_sh = ImageTranslateFFT(Im, dx, dy)
% Im_sh = ImageTranslateFFT(Im, dx, dy)
%
% dx, dy are pixels

[Ny, Nx] = size(Im);

[~, ~, X, Y] = CreateGrid([Ny, Nx]);

% fft and multiply by linear phase
F = fft2(Im);
F_sh = F .* fftshift( exp(1i*2*pi*(X*dx/Nx + Y*dy/Ny)) );

Im_sh = ifft2(F_sh);

if isreal(Im),
    Im_sh = real(Im_sh);
end

