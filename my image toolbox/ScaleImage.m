function Imout = ScaleImage(Img, ascale, varargin)
% Imout = ScaleImage(Img, ascale, varargin)
%
% use fft then DFT to scale an image around the mid-point of the array
%
% assumes Img is real

bReal = CheckOption('isreal', true, varargin{:});

[x, y] = CreateGrid(Img);
[nr, nc] = size(Img);

% FFT
IMG = fftshift(fft2(ifftshift(Img)));

% DFT Matrices using outer products and scale
PX = exp(1j*(2*pi/ascale)*(x*x')/nc);
PY = exp(1j*(2*pi/ascale)*(y*y')/nr);
Imout = (PX*IMG)*PY/(nr*nc);

if bReal
    Imout = real(Imout);
end

%figure, imageschcit(abs(Imout - Img)), colorbar
