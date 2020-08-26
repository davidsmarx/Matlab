import matlab.io.*
if exist(out_fn, 'file'), delete(out_fn); end
fptr = fits.createFile(out_fn);
fits.createImg(fptr,'double_img',size(pha_init));
fits.writeImg(fptr, pha_init);
fits.writeKey(fptr, 'pupdiam', 2*sOptionsscale.Rnorm, 'pupil diameter (pix)');
fits.createImg(fptr,'double_img',size(mask_init))
fits.writeImg(fptr, double(mask_init));
fits.closeFile(fptr);
