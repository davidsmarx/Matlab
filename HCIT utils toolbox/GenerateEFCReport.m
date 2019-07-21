function sOut = GenerateEFCReport(runnum, listItnum, csDisplayFun, varargin)
% sOut = GenerateEFCReport(runnum, listItnum, sDisplayFun)
%
% listItnum can be array of itnums, or array of CRunData objects
% csDisplayFun is a cell array of methods, e.g.{'DisplayAllInt','DisplayCEfields'}
%
% create figures S.(sDisplayFun) for each iteration in listItnum
% and copy to a PowerPoint presentation
% 
% 

%mlock   
persistent Sppt;

% options
ppt_fn = CheckOption('pptfn', '', varargin{:});
    
% check if PowerPoint Presentation already exists and still there
try
    Sppt.Presentation.Slides,
catch
    clear Sppt; Sppt = [];
end

% open PowerPoint if necessary
if isempty(Sppt) && ispc,
    Sppt = Cppt(ppt_fn);
end

% get the CRunData objects
N = length(listItnum);
if isnumeric(listItnum),
    for ii = 1:N
        fprintf('reading itnum %d\n', listItnum(ii));
        S(ii) = CRunData(runnum, listItnum(ii));
    end
elseif isa(listItnum, 'CRunData')
    S = listItnum;
else
    error(['listItnum type error: ' class(listItnum)]);
end
  
for iplot = 1:length(csDisplayFun),
    listHfig{iplot} = CreaetePlots(S, csDisplayFun{iplot}, Sppt, varargin{:});
end

    if nargout >= 1,
        sOut = struct(...
            'listS', S ...
            ,'listHfig', listHfig ...
            ,'Sppt', Sppt ...
            );
    end

end % main

function hfig = CreaetePlots(S, sDisplayFun, Sppt, varargin)
    % create the plots
    % some plots are differential

    save_pn = CheckOption('save_pn', ['./' sDisplayFun '/'], varargin{:}); % is ~ispc
    figheight = CheckOption('figheight', 700, varargin{:}); % for ppt display

    % create path to put plots, if necessary
    if ~ispc && ~exist(save_pn)
        mkdir(save_pn)
    end
    
    % list of CRunData methods where the first argument is a reference
    % iteration:
    listDiff = {
        'DisplayDEfields'
        'DisplayDMv'
        'DisplayCEfields'
        };
    
    N = length(S);

    if any(strcmp(sDisplayFun, listDiff)),
        for ii = 1:N-1,
            hfig(ii) = S(ii+1).(sDisplayFun)(S(ii), varargin{:});
            figscale = CalcFigscale(hfig(ii), figheight);
            set(hfig(ii), 'Position', figscale*get(hfig(ii),'position'));
            if ispc,
                htmp = Sppt.CopyFigNewSlide(hfig(ii));
                %set(htmp,'Height',figheight);
            else
                saveas(hfig(ii), [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
            end
        end
        
    else, % one call per iteration
        for ii = 1:N,
            hfig(ii) = S(ii).(sDisplayFun)(varargin{:});
            figscale = CalcFigscale(hfig(ii), figheight);
            set(hfig(ii), 'Position', figscale*get(hfig(ii),'position'));
            if ispc,
                htmp = Sppt.CopyFigNewSlide(hfig(ii));
                %set(htmp,'Height',figheight);
            else
                saveas(hfig(ii), [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
            end
        end
    end

end % CreatePlots

function figscale = CalcFigscale(hfig, figheight)
    
    pos = get(hfig,'Position');
    ysize = pos(end);
    figscale = figheight/ysize;

end % CalcFigscale