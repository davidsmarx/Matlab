function haxout = imagescampphase(A,x,y,varargin)
% hax = imagescampphase(A,x,y,options)
%
% x,y can be [], scalar dx(y), or vector x, y
%
% options:
%   'bLog' or 'scale', (false) or 'true' (log10) or 'log10' or 'db'
%   'title', titlestr
%   ydir = CheckOption('ydir', 'reverse', varargin{:}); or 'normal'
%   xlabelstr = CheckOption('xlabel', [], varargin{:});
%   ylabelstr = CheckOption('ylabel', [], varargin{:});
%   xylim = CheckOption('xylim', [], varargin{:});
%   phasescale = CheckOption('phasescale', 'pi rad', varargin{:}); % or 'rad', 'deg'

[ny, nx] = size(A);

if nargin < 2 || isempty(x),
    x = (-nx/2:nx/2-1)';
end
if nargin < 3 || isempty(y),
    y = (-ny/2:ny/2-1)';
end
if isscalar(x),
    x = x*(-nx/2:nx/2-1)';
end
if isscalar(y),    
    y = y*(-ny/2:ny/2-1)';
end

% options
titlestr = CheckOption('title','',varargin{:});
bLog = CheckOption('bLog',false,varargin{:});
if ~bLog, bLog = CheckOption('scale',false,varargin{:}); end % alternate option keyword
ydir = CheckOption('ydir', 'reverse', varargin{:});
xlabelstr = CheckOption('xlabel', [], varargin{:});
ylabelstr = CheckOption('ylabel', [], varargin{:});
xylim = CheckOption('xylim', [], varargin{:});
phasescale = CheckOption('phasescale', 'pi rad', varargin{:}); % or 'rad', 'deg'

pos = get(gcf,'position');
pos(3) = 1.75*pos(3);
set(gcf,'position',pos)

hax(1) = subplot(1,2,1);
if islogical(bLog) && bLog,
    bLog = 'log10';
end
switch lower(bLog)
    case {'log','log10'},
        Aplot = real(log10(abs(A).^2));
        strColorbarAmp = 'log_{10} Intensity';
    case 'db',
        Aplot = real(20*log10(abs(A)));
        strColorbarAmp = '|A|^2 (dB)';
    otherwise
        % 'linear' is default
        Aplot = abs(A);
        strColorbarAmp = '|A|';
end

% whether to plot rad or pi rad
switch lower(phasescale),
    case 'pi rad'
        phasenorm = pi;
        sCtitle = '\pi rad';
    case 'rad'
        phasenorm = 1;
        sCtitle = 'rad';
    case 'deg'
        phasenorm = pi/180;
        sCtitle = 'deg';
    otherwise
        error(['unknown phasescale option: ' phasescale]);
end

imagesc(x,y,Aplot), axis image,...
    title([pwd2titlestr(titlestr) ', Amp']), SetColorbar(strColorbarAmp);
if ~isempty(xlabelstr), xlabel(xlabelstr), end
if ~isempty(ylabelstr), ylabel(ylabelstr), end


hax(2) = subplot(1,2,2);
imagesc(x,y,angle(A)./phasenorm), axis image, title([pwd2titlestr(titlestr) ', Phase']),...
    colorbartitle(sCtitle) %SetColorbar('\pi radians');
if ~isempty(xlabelstr), xlabel(xlabelstr), end
if ~isempty(ylabelstr), ylabel(ylabelstr), end

set(hax,'ydir',ydir);

if ~isempty(xylim),
    set(hax,'xlim',xylim,'ylim',xylim)
end

if nargout > 0,
    haxout = hax;
end

return

function htit = SetColorbar(titlestr)

h_cb = colorbar;
htit = get(h_cb,'title');
set(htit,'string',titlestr);

return