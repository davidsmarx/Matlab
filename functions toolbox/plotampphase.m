function [hamp, hpha] = plotampphase(x,y,varargin)
% [hamp, hph] = plotampphase(x,y,optionslist)
% Supported options are:
%   'dB'
%   'title', followed by a string
%   'xlabel', followed by a string
%   'legend', followed Npl strings where Npl is the # of columns of y
%   'xlim', 'ylim', followed by two element vector of the limits
%
% return:
%   hamp = handle to amplitude axes
%   hph  = handle to phase axes

if nargin == 0, error('usage: han = plotampphase(x,y,options)'); end

if nargin == 1, y = x; x = [0:length(y)-1]; end

%set(gcf,'position',[360 1 560 560]),
figure_mxn(gcf, 2, 1);

[nab, npl] = size(y);
% initialize options
yplot = abs(y);
ylab  = 'Magnitude [linear]';
ptitle = [];
xlab = [];
leg = 'hide';
xlimprop = 'xlimmode'; xlim = 'auto';
ylimprop = 'ylimmode'; ylim = 'auto';
if ~isempty(varargin),
	
	if strmatch('dB',varargin), 
		yplot = 2*decibel(y);
		ylab  = 'Magnitude [dB]';
	end

	itit = strmatch('title',varargin);
	if itit, ptitle = varargin{itit+1}; end
	
	ixla = strmatch('xlabel',varargin);
	if ixla, xlab = varargin{ixla+1}; end
	
	ileg = strmatch('legend',varargin);
	if ileg, 
		leg = cell(1,npl); 
		[leg{:}] = deal(varargin{ileg+1:ileg+npl}); 
	end

    ixlim = strmatch('xlim',varargin);
    if ixlim,
        xlim = str2num(varargin{ixlim+1});
        xlimprop = 'xlim';
    end
    
    iylim = strmatch('ylim',varargin);
    if iylim,
        ylim = str2num(varargin{iylim+1});
        ylimprop = 'ylim';
    end

end

h_amp = subplot(2,1,1);
plot(x,yplot),
title(ptitle),
ylabel(ylab),
xlabel(xlab),
legend(leg),
grid
set(gca,ylimprop,ylim), set(gca,xlimprop,xlim)

h_pha = subplot(2,1,2);
plot(x,angle(y)/pi),
ylabel('[\pi rad]'),
xlabel(xlab),
grid,
set(gca,ylimprop,ylim), set(gca,xlimprop,xlim)

if nargout > 0,
    hamp = h_amp;
    hpha = h_pha;
end
    
return
