function sOut = GenerateCGSReport(listGsnum, bn, varargin)
% sOut = GenerateCGSReport(listGsnum, bn, cDisplayFun)
%
% listGsnum is array of gsnums
% bn is same as bn in CGS()
% varargin = additional cdDisplayFun cell arrays as you want
%
% create figures S.(sDisplayFun) for each iteration in listGsnum
% and copy to a PowerPoint presentation (Windows) or save .png to results
% folder
% 
% methods:
% DisplayGS
% DisplayGSrefGS
% ZernZernRemapFit
% DisplayAllPlanes

more off

% mlock   
% persistent Sppt;

% options
ppt_fn = CheckOption('pptfn', '', varargin{:});
%sOptin = CheckOption('sOptin', [], varargin{:}); % passed to CGS
Sppt = CheckOption('Sppt', [], varargin{:});

% % check if PowerPoint Presentation already exists and still there
% try
%     Sppt.Presentation.Slides,
% catch
%     clear Sppt; Sppt = [];
% end

% open PowerPoint if necessary, and plots are requested
if isempty(Sppt) && ispc && ~isempty(varargin),
    Sppt = Cppt(ppt_fn);
end

% get the CGS objects
N = length(listGsnum);
if isnumeric(listGsnum),
    for ii = 1:N
        fprintf('reading gsnum %d\n', listGsnum(ii));
        S(ii) = CGS(listGsnum(ii), bn); %, sOptin);
    end
elseif isa(listItnum, 'CGS')
    S = listItnum;
else
    error(['listItnum type error: ' class(listItnum)]);
end

% call the plotting methods
if length(varargin) == 0, listHfig = []; end
for iplot = 1:length(varargin),
    if iscell(varargin{iplot}),
        listHfig(iplot) = CreatePlots(S, varargin{iplot}{1}, Sppt, varargin{iplot}{2:end});
    end
end

% plot graphs of metrics v itnum

more on

sOut = S;

end % main

function [hfig] = CreatePlots(S, sDisplayFun, Sppt, varargin)
    % create the plots
    % some plots are differential
    % some plots we also plot metrics v itnum

    save_pn = CheckOption('save_pn', ['./' sDisplayFun '/'], varargin{:}); % is ~ispc
    figheight = CheckOption('figheight', 700, varargin{:}); % for ppt display
    trialname = CheckOption('trialname', '', varargin{:});

    % create path to put plots, if necessary
    if ~ispc && ~exist(save_pn)
        mkdir(save_pn)
    end
    
    %     % list of CRunData methods where the first argument is a reference
    %     % iteration:
    %     listRef = {
    %         'DisplayGSrefGS'
    %         };
    
    N = length(S);

    for ii = 1:N,
        S(ii).(sDisplayFun)(varargin{:}); %,'hfig',hfig);
        hfig = gcf;
        figscale = CalcFigscale(hfig, figheight);
        set(hfig, 'Position', figscale*get(hfig,'position'));
        if ispc,
            htmp = Sppt.CopyFigNewSlide(hfig);
            %set(htmp,'Height',figheight);
        else
            saveas(hfig, [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
        end
    end
    

end % CreatePlots

function figscale = CalcFigscale(hfig, figheight)
    
    pos = get(hfig,'Position');
    ysize = pos(end);
    figscale = figheight/ysize;

end % CalcFigscale

