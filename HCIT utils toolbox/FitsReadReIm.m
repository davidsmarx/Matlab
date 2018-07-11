function Im = FitsReadReIm(bn)
% Im = FitsReadReIm(bn)

Im = fitsread([bn 'real.fits']) + 1i*fitsread([bn 'imag.fits']);

