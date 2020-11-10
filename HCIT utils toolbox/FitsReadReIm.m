function Im = FitsReadReIm(bn)
% Im = FitsReadReIm(bn)
% 
% either
%    fitsread(bn) + 1j*fitsread(bn,'image');
% or
%    fitsread([bn 'real.fits']) + 1i*fitsread([bn 'imag.fits']);

if exist(bn,'file'),
    Im = fitsread(bn) + 1j*fitsread(bn,'image');
    
else,
    Im = fitsread([bn 'real.fits']) + 1i*fitsread([bn 'imag.fits']);
    
end

