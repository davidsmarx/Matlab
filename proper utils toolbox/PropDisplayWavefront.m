function [hfig, hax] = PropDisplayWavefront(wavefront, varargin)
% [hfig, hax] = PropDisplayWavefront(wavefront, varargin)
%
% wavefront is a Proper wavefront

U = CConstants;

ampclim = CheckOption('ampclim', [], varargin{:});
phaclim = CheckOption('phaclim', [-1 1], varargin{:});
hfig = CheckOption('hfig', [], varargin{:});
%hax = CheckOption('hax', [], varargin{:});
strtitle = CheckOption('title', ' ', varargin{:});

[x, y] = CreateGrid(wavefront.wf, wavefront.dx);

if isempty(hfig), hfig = figure; else, figure(hfig); end

E = fftshift(wavefront.wf);

%Img = abs(E).^2;
%imageschcit(x/U.MM, y/U.MM, real(log10(Img))), colorbartitle('log_{10} Intensity')
hax = imagescampphase(E, x/U.MM, y/U.MM, 'bLog', true, 'ydir', 'normal','title',strtitle);
xlabel('X (mm)'), ylabel('Y (mm)')
if ~isempty(ampclim), set(hax(1),'clim',ampclim); end
if ~isempty(phaclim), set(hax(2),'clim',phaclim); end


