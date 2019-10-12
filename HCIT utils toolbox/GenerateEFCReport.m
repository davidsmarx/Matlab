function sOut = GenerateEFCReport(runnum, listItnum, varargin)
% sOut = GenerateEFCReport(runnum, listItnum, cDisplayFun)
%
% listItnum can be array of itnums, or array of CRunData objects
% csDisplayFun is a cell array: {method, varargin options (e.g. 'clim', clim)}
% varargin = additional cdDisplayFun cell arrays as you want
%
% create figures S.(sDisplayFun) for each iteration in listItnum
% and copy to a PowerPoint presentation (Windows) or save .png to results
% folder
% 
% methods:
% DisplayImCubeImage
% DisplayImCubeUnProb
% DisplayImCubeContrast
% DisplayImCubeSigProb
% DisplayIncInt
% DisplayProbeAmp
% DisplayProbeCube
% DisplayCohInt
% DisplayAllInt
% DisplayRadialIntensity
% DisplayIncCohInt
% DisplayEfields
% DisplayDEfields
% DisplayCEfields
% DisplayDMv
% DisplayDMvProbe

more off

%mlock   
persistent Sppt;

% options
ppt_fn = CheckOption('pptfn', '', varargin{:});
sOptin = CheckOption('sOptin', [], varargin{:});
    
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
        S(ii) = CRunData(runnum, listItnum(ii), sOptin);
    end
elseif isa(listItnum, 'CRunData')
    S = listItnum;
else
    error(['listItnum type error: ' class(listItnum)]);
end

listHfig = {};
for iplot = 1:length(varargin),
    if iscell(varargin{iplot}),
        listHfig{end+1} = CreaetePlots(S, varargin{iplot}{1}, Sppt, varargin{iplot}{2:end});
    end
end

% plot graphs of metrics v itnum
[hfig, haxprobeh, probeh] = PlotProbeh(S);
[hfig, haxrmsddmv, rmsdDMv] = PlotRMSdDMv(S);

    if nargout >= 1,
        sOut = struct(...
            'listS', S ...
            ,'listHfig', listHfig ...
            ,'Sppt', Sppt ...
            ,'probeh', probeh ...
            ,'haxprobeh', haxprobeh ...
            ,'haxrmsddmv', haxrmsddmv ...
            ,'rmsdDMv', rmsdDMv ...
            );
    end

more on

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
        };
    
    N = length(S);

    hfig = [];
    switch sDisplayFun,
        case listDiff
            %if any(strcmp(sDisplayFun, listDiff)),
            for ii = 1:N-1,
                hfig = S(ii+1).(sDisplayFun)(S(ii), varargin{:},'hfig',hfig);
                figscale = CalcFigscale(hfig, figheight);
                set(hfig, 'Position', figscale*get(hfig,'position'));
                if ispc,
                    htmp = Sppt.CopyFigNewSlide(hfig);
                    %set(htmp,'Height',figheight);
                else
                    saveas(hfig, [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
                end
            end % for ii iter
        
        case 'DisplayCEfields'

            for ii = 1:N-1,
                [hfig, hax, sCtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:},'hfig',hfig);
                sCmetric(ii) = sCtmp;
                
                figscale = CalcFigscale(hfig, figheight);
                set(hfig, 'Position', figscale*get(hfig,'position'));
                if ispc,
                    htmp = Sppt.CopyFigNewSlide(hfig);
                    %set(htmp,'Height',figheight);
                else
                    saveas(hfig, [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
                end
            end % for ii iter

            figure, plotampphase([S(2:N).iter], [sCmetric.CC],...
                'xlabel','Iteration #','title','\DeltaE Testbed Model Correlation (CC)');
            
        otherwise, % one call per iteration
            for ii = 1:N,
                hfig = S(ii).(sDisplayFun)(varargin{:},'hfig',hfig);
                figscale = CalcFigscale(hfig, figheight);
                set(hfig, 'Position', figscale*get(hfig,'position'));
                if ispc,
                    htmp = Sppt.CopyFigNewSlide(hfig);
                    %set(htmp,'Height',figheight);
                else
                    saveas(hfig, [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
                end
            end
    end

end % CreatePlots

function figscale = CalcFigscale(hfig, figheight)
    
    pos = get(hfig,'Position');
    ysize = pos(end);
    figscale = figheight/ysize;

end % CalcFigscale

function [hfig, hax, probeh] = PlotProbeh(S)

     [itnum, probeh] = deal(zeros(size(S)));
     for ii = 1:length(S)
         itnum(ii)  = S(ii).iter;
         probeh(ii) = FitsGetKeywordVal(S(ii).ImKeys,'PROBEH');
     end

     hfig = figure;
     semilogy(itnum, probeh, '-o'), grid
     xlabel('Iteration #')
     ylabel('probeh')
     hax = gca;
     
end % PlotProbeh

function [hfig, hax, rmsdDMv] = PlotRMSdDMv(listS, varargin)

    for ii = 1:length(listS),
        if isempty(listS(ii).DMvCube)
            listS(ii).ReadDMvCube;
        end
    end
    
    Ndm = length(listS(1).DMvCube);
    
    for ii = 1:length(listS),
        for idm = 1:Ndm
            DMvtmp{idm} = squeeze(listS(ii).DMvCube{idm}(:,:,1));
        end
        listDMv{ii} = [DMvtmp{:}];
    end
    
    % rms difference
    rmsdDMv = zeros(length(listS)-1,1);
    for ii = 1:length(listS)-1
        iuse = abs(listDMv{ii+1}) > 0 & abs(listDMv{ii}) > 0;
        rmsdDMv(ii) = rms(listDMv{ii+1}(iuse) - listDMv{ii}(iuse));
    end
    
    hfig = figure;
    plot([listS(2:end).iter].', rmsdDMv), grid
    xlabel('Iteration #')
    ylabel('rms \Delta Vmu')
    hax = gca;
    
end % PlotRMSdDMv