function [rplot, Ar, ha, hl] = plotradial(xi, yi, A, varargin)
% [hfig, ha, hl, rplot, IntRad] = plotradial(x, y, A, varargin)
%
% generic routine for plotting azimuthal average of A vs radius
%
% % A is an array or cell array of arrays
%
% return [rplot, Ar, ha, hl]
%
%             Nr = CheckOption('nr', ceil(min([128 length(R)/4])), varargin{:}); % # of radial sample pts
%             dispRadlim = CheckOption('dispradlim', [0 S.XYlimDefault], varargin{:});
%             drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
%             bMaskUse = CheckOption('bMask', S.bMask, varargin{:});
%             strYlabel = CheckOption('ylabel', 'Average Normalized Intensity', varargin{:});
%             plotRequired = CheckOption('plotrequired', [], varargin{:}); % [r(:) contrast(:)]
%             iplot = CheckOption('iplot', 1:length(ImCube), varargin{:});
%             legstr = CheckOption('legstr', [], varargin{:});
%             strTitle = CheckOption('title', ['Iter #' num2str(S.iter)], varargin{:});
%             bPlotMean = CheckOption('plotmean', true, varargin{:});
%             haxuse = CheckOption('hax', [], varargin{:});

% coordinates and checks
[nr, nc] = size(A);
x = xi; y = yi; % default
if isempty(xi) || isempty(yi),
    [x, y] = CreateGrid(A);
end
if isscalar(xi), x = xi + (0:(nc-1)).'; end
if isscalar(yi), y = yi + (0:(nr-1)).'; end
if length(xi) ~= nc || length(yi) ~= nr
    error('size mismatch');
end
[X, Y] = meshgrid(x, y); R = hypot(X, Y);

% options
Nr = CheckOption('nr', ceil(min([128 length(R)/4])), varargin{:}); % # of radial sample pts
dispRadlim = CheckOption('dispradlim', [0 max(abs([X(:); Y(:)]))], varargin{:});
scale = CheckOption('scale', 'linear', varargin{:}); % or 'log'
drawRadii = CheckOption('drawradii', [], varargin{:}); % draw vertical lines
bMaskUse = CheckOption('bMask', 1, varargin{:});
strXlabel = CheckOption('xlabel', 'Radius', varargin{:});
strYlabel = CheckOption('ylabel', 'Aximuthal Average', varargin{:});
plotRequired = CheckOption('plotrequired', [], varargin{:}); % [r(:) contrast(:)]
legstr = CheckOption('legstr', [], varargin{:});
strTitle = CheckOption('title', [], varargin{:});
haxuse = CheckOption('hax', [], varargin{:});

%
re = linspace(dispRadlim(1), dispRadlim(2), Nr+1)';
Ar = zeros(Nr,1);
for ir = 1:Nr,
    % including S.bMask applies theta (bowtie) limits
    Ar(ir) = mean(A(R > re(ir) & R <= re(ir+1) & bMaskUse));
end % for ir
rplot = mean([re(1:end-1) re(2:end)],2); % radii midway between edges

%if isempty(legstr), legstr = legstrwv; end

if ~isempty(haxuse),
    axes(haxuse);
end

switch scale
    case 'log'
        hl = semilogy(rplot, Ar);
    case 'linear'
        hl = plot(rplot, Ar);
    otherwise
        error(['unknown scale ' scale]);
end

ha = gca;
hold on


grid on

if ~isempty(drawRadii),
    ylim = get(gca,'ylim');
    hold on
    for irad = 1:length(drawRadii),
        plot(drawRadii(irad)*[1 1], ylim, '--r')
    end
    hold off
end

xlabel(strXlabel)
ylabel(strYlabel)
%hleg = legend(legstr{:}) %, 'location','north');
title(strTitle)


end % DisplayRadialPlot
