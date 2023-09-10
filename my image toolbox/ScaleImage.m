function Imout = ScaleImage(Img, ascale, varargin)
% Imout = ScaleImage(Img, ascale, varargin)
%
% use fft then DFT to scale an image around the mid-point of the array
% if scale < 1, then output image is aliased. Default behavior is to set
% the aliased part of the output image to zero
%
% assumes Img is real
% options:
%   'isreal', (true) or false

bReal = CheckOption('isreal', true, varargin{:});

[x, y, X, Y] = CreateGrid(Img);
[nr, nc] = size(Img);

% FFT
IMG = fftshift(fft2(ifftshift(Img)));

% DFT Matrices using outer products and scale
PX = exp(1j*(2*pi/ascale)*(x*x')/nc);
PY = exp(1j*(2*pi/ascale)*(y*y')/nr);
bNotAlias = abs(X) < 0.5*nc*ascale & abs(Y) < 0.5*nr*ascale;

Imout = bNotAlias .* ((PX*IMG)*PY/(nr*nc));


if bReal
    Imout = real(Imout);
end

%figure, imageschcit(abs(Imout - Img)), colorbar
