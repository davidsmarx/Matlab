function h = imageWavefrontUNF(unfname)
% h = imageWavefrontUNF(unfname)

MM = 1e-3;

if nargin == 0,
    [filename, pathname, fi] = uigetfile('*.unf','Select a UNF file');
    unfname = [pathname filename];
end

[ A, dx, dy, lambda, curv, xp, yp ] =...
    ReadWavefrontUNF( unfname );

h = figure;
imagesc(xp/MM,yp/MM,abs(A)), axis image,
title(pwd2titlestr(unfname,0));