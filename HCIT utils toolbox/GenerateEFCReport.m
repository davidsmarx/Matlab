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

% mlock   
% persistent Sppt;

% options
ppt_fn = CheckOption('pptfn', '', varargin{:});
sOptin = CheckOption('sOptin', [], varargin{:}); % passed to CRunData
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

% call the plotting methods
if length(varargin) == 0, listHfig = []; end
for iplot = 1:length(varargin),
    if iscell(varargin{iplot}),
        listHfig(iplot) = CreaetePlots(S, varargin{iplot}{1}, Sppt, varargin{iplot}{2:end});
    end
end

% plot graphs of metrics v itnum
[hfig, haxprobeh, probeh] = PlotProbeh(S);
[hfig, haxrmsddmv, rmsdDMv] = PlotRMSdDMv(S);
[hfig, haxtexp, texp] = PlotTexp(S);

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

function [hfig, hax, sMetrics] = CreaetePlots(S, sDisplayFun, Sppt, varargin)
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
                [hfig, hax, sMtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:},'hfig',hfig);
                sMetrics(ii) = sMtmp;
                
                figscale = CalcFigscale(hfig, figheight);
                set(hfig, 'Position', figscale*get(hfig,'position'));
                if ispc,
                    htmp = Sppt.CopyFigNewSlide(hfig);
                    %set(htmp,'Height',figheight);
                else
                    saveas(hfig, [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
                end
            end % for ii iter
        
            if strcmp(sMtmp(1).type, 'dEfields'),
                nw = S(1).NofW; % for convenience
                figure, 
                hh = semilogy([S(2:N).iter], [sMetrics.rmsdE_t].^2, '-', ...
                    [S(2:N).iter], mean([sMetrics.rmsdE_t].^2, 1), '-', ...
                    [S(2:N).iter], [sMetrics.rmsdE_m].^2, '--', ...
                    [S(2:N).iter], mean([sMetrics.rmsdE_m].^2, 1), '--');
                
                grid on
                set(hh(nw+1),'linewidth',2)
                set(hh(end), 'linewidth',2)
                %set(gca,'ylim',get(hax(1),'clim'))
                %set(gca,'ylim',[1e-9 1e-6])
                wvstrtmp = num2str(S(1).NKTcenter(1:nw)'/S(1).NM);
                tbwvstrtmp = char(strcat('TB', {' '}, wvstrtmp, 'nm'));
                mowvstrtmp = char(strcat('Model', {' '}, wvstrtmp, 'nm'));
                legend(char(tbwvstrtmp, 'TB-mean', mowvstrtmp, 'Model-Mean'))
                xlabel('Iteration #')
                ylabel('mean |\DeltaE|^2')
                title(trialname, 'fontsize', 14)
            end
            
        case 'DisplayCEfields'

            for ii = 1:N-1,
                [hfig, hax, sCtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:},'hfig',hfig);
                sMetrics(ii) = sCtmp;
                
                figscale = CalcFigscale(hfig, figheight);
                set(hfig, 'Position', figscale*get(hfig,'position'));
                if ispc,
                    htmp = Sppt.CopyFigNewSlide(hfig);
                    %set(htmp,'Height',figheight);
                else
                    saveas(hfig, [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
                end
            end % for ii iter

            figure, plotampphase([S(2:N).iter], [sMetrics.CC],...
                'xlabel','Iteration #','title',[trialname ', \DeltaE Testbed Model Correlation (CC)']);
            
        otherwise, % one call per iteration
            sMetrics = struct;
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

function [hfig, hax, texp] = PlotTexp(S)

    [itnum, texp] = deal(zeros(size(S)));
    for ii = 1:length(S)
        itnum(ii) = S(ii).iter;
        texp(ii) = FitsGetKeywordVal(S(ii).ImKeys, 'texp');
    end
    
    hfig = figure;
    plot(itnum, texp, 'o'), grid
    xlabel('Iteration #')
    ylabel('T_{exp} (s)')
    hax = gca;
    
end % PlotTexp

function [hfig, hax, rmsdDMv] = PlotRMSdDMv(listS, varargin)
    % [hfig, hax, rmsdDMv] = PlotRMSdDMv(listS, varargin)

    % read all the DMv cubes
    for ii = 1:length(listS),
        if isempty(listS(ii).DMvCube)
            listS(ii).ReadDMvCube;
        end
    end
    
    Ndm = length(listS(1).DMvCube);
    
    % select the unprobed DMv from each cube
    for ii = 1:length(listS),
        for idm = 1:Ndm
            DMvtmp{idm} = squeeze(listS(ii).DMvCube{idm}(:,:,1));
        end
        listDMv(ii,:) = DMvtmp; % listDMv is cell array
    end
    
    % rms difference
    rmsdDMv = zeros(length(listS)-1,Ndm);
    for ii = 1:length(listS)-1
        for idm = 1:Ndm,
            iuse = abs(listDMv{ii+1,idm}) > 0 & abs(listDMv{ii,idm}) > 0;
            rmsdDMv(ii,idm) = rms(listDMv{ii+1,idm}(iuse) - listDMv{ii,idm}(iuse));
        end
    end
    
    hfig = figure;
    hh = plot([listS(2:end).iter].', rmsdDMv, '-', [listS(2:end).iter].', mean(rmsdDMv,2), '--');
    set(hh,'LineWidth', 1.0)
    grid on
    xlabel('Iteration #')
    ylabel('rms \Delta Vmu')
    for idm = 1:Ndm,
        legstr{idm} = ['DM ' num2str(idm)];
    end
    legend(legstr{:}, 'Mean')
    hax = gca;
    
end % PlotRMSdDMv