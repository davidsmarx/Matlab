function [hfig, hax] = PropDisplayWavefront(wavefront, varargin)
% [hfig, hax] = PropDisplayWavefront(wavefront, varargin)
%
% wavefront is a Proper wavefront

U = CConstants;

[x, y] = CreateGrid(wavefront.wf, wavefront.dx);

hfig = figure;

E = fftshift(wavefront.wf);
Img = abs(E).^2;

imageschcit(x/U.MM, y/U.MM, real(log10(Img))), colorbartitle('log_{10} Intensity')
set(gca,'clim',[-9 0])

hax = gca;
