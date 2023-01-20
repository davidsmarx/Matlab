% NAME:
%     fft_rotate
%
% PURPOSE:
%     Shift and rotate an image in one step via FFTs
%
% EXPLANATION:
%     fft_rotate first rotates an image by ANGLE degrees about the image center
%     and then shifts it by xshift,yshift pixels.  The image must be square and
%     should be a factor of two in dimension and zero-padded, since FFTs are used.
%
% CALLING SEQUENCE:
%     result = fft_rotate( image, angle [, xshift, yshift] )
%
% INPUT PARAMETERS:
%     image : 2-d square image array (may be complex-valued)
%     angle : Angle in degrees counterclockwise to rotate image
%
% OPTIONAL INPUT PARAMETERS:
%     xshift, yshift : Number of pixels to shift rotated image (may be subpixel)
%
% RESULT:
%     Returns complex-valued rotated and shifted image
%
% AUTHOR:
%     John Krist, based on a method presented by Larkin et al. in Optics Communications (1997)

function g = fft_rotate( image, t_deg, xshift, yshift )

if nargin < 4; yshift = 0; end
if nargin < 3; xshift = 0; end

[ny, dim] = size(image);	% assume square image

t = -t_deg / 180. * pi;
a = tan(t/2.);
b = -sin(t);

x  = [-floor(dim / 2 ) : floor((dim - 1) / 2)];
[x, y] = meshgrid(x, x);
u = circshift( x / dim, -dim/2, 2 );

gx =  ifft(exp(-2*pi*1i*u.*a.*y) .* fft(image, dim, 2), dim, 2);
gyx = ifft(exp(-2*pi*1i*transpose(u).*(b*x+yshift)) .* fft(gx, dim, 1), dim, 1);
g =   ifft(exp(-2*pi*1i*u.*(a*(y-yshift)+xshift)) .* fft(gyx, dim, 2), dim, 2);

end
