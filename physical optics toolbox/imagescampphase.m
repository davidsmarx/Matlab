function haxout = imagescampphase(A,x,y,varargin)
% hax = imagescampphase(A,x,y,options)
%
% options:
%   'bLog', true or (false)
%   'titlestr', titlestr

[ny, nx] = size(A);

if nargin < 2 || isempty(x),
    x = (-nx/2:nx/2-1)';
end
if nargin < 3 || isempty(y),
    y = (-ny/2:ny/2-1)';
end

% options
titlestr = CheckOption('titlestr','',varargin{:});
bLog = CheckOption('bLog',false,varargin{:});
ydir = CheckOption('ydir', 'reverse', varargin{:});
xlabelstr = CheckOption('xlabel', [], varargin{:});
ylabelstr = CheckOption('ylabel', [], varargin{:});

pos = get(gcf,'position');
pos(3) = 1.75*pos(3);
set(gcf,'position',pos)

hax(1) = subplot(1,2,1);
if bLog,
    Aplot = real(20*log10(abs(A)));
    strColorbarAmp = '|A|^2 (dB)';
else,
    Aplot = abs(A);
    strColorbarAmp = '|A|';
end

imagesc(x,y,Aplot), axis image,...
    title([pwd2titlestr(titlestr) ', Amp']), SetColorbar(strColorbarAmp);
if ~isempty(xlabelstr), xlabel(xlabelstr), end
if ~isempty(ylabelstr), ylabel(ylabelstr), end


hax(2) = subplot(1,2,2);
imagesc(x,y,angle(A)./pi), axis image, title([pwd2titlestr(titlestr) ', Phase']),...
    colorbartitle('\pi rad') %SetColorbar('\pi radians');
if ~isempty(xlabelstr), xlabel(xlabelstr), end
if ~isempty(ylabelstr), ylabel(ylabelstr), end

set(hax,'ydir',ydir);

if nargout > 0,
    haxout = hax;
end

return

function htit = SetColorbar(titlestr)

h_cb = colorbar;
htit = get(h_cb,'title');
set(htit,'string',titlestr);

return