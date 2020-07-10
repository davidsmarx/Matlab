function [hfig_out, hax_out] = PropDisplayWavefront(wavefront, varargin)
% [hfig, hax] = PropDisplayWavefront(wavefront, varargin)
%
% wavefront is a Proper wavefront
%
% Options:
%   scale = CheckOption('scale', 'log', varargin{:}); % or linear
%   ampclim = CheckOption('ampclim', [], varargin{:});
%   phaclim = CheckOption('phaclim', [-1 1], varargin{:});
%   hfig = CheckOption('hfig', [], varargin{:});
%   %hax = CheckOption('hax', [], varargin{:});
%   strtitle = CheckOption('title', ' ', varargin{:});

U = CConstants;

scale = CheckOption('scale', 'log', varargin{:}); % or linear
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
hax = imagescampphase(E, x/U.MM, y/U.MM, 'bLog', strcmpi(scale,'log'), 'ydir', 'normal','title',strtitle, ...
    'xlabel', 'X (mm)', 'ylabel', 'Y (mm)');
if ~isempty(ampclim), set(hax(1),'clim',ampclim); end
if ~isempty(phaclim), set(hax(2),'clim',phaclim); end

if nargout > 0,
    hfig_out = hfig;
    hax_out = hax;
end

