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
    %if isscalar(listItnum)
    %ii = 1;
    %while true,
    for ii = 1:N;
        fprintf('reading itnum %d\n', listItnum(ii));
        %try
            S(ii) = CRunData(runnum, listItnum(ii), sOptin);
        %catch
        %    warning(['iter# ' num2str(listItnum(ii)) ' not found']);
        %end
        
    end
elseif isa(listItnum, 'CRunData')
    S = listItnum;
else
    error(['listItnum type error: ' class(listItnum)]);
end

% plot graphs of metrics v itnum on the first slide
slide = Sppt.NewSlide(1);
hfig = figure_mxn(2,2);
hax(1,1) = subplot(2,2,1); [~, ~, ~, itnum_min] = PlotNormIntensity(S, 'hfig', hfig, 'hax', hax(1,1));
hax(1,2) = subplot(2,2,2); PlotBeta(S, 'hfig', hfig, 'hax', hax(1,2));
hax(2,1) = subplot(2,2,3); probeh = PlotProbeh(S, 'hfig', hfig, 'hax', hax(2,1));
hax(2,2) = subplot(2,2,4); rmsdDMv = PlotRMSdDMv(S, 'hfig', hfig, 'hax', hax(2,2));
hPic = Sppt.CopyFigSlide(slide, hfig);

% call the plotting methods
listHfig = {};
for iplot = 1:length(varargin),
    if iscell(varargin{iplot}),
        listHfig{end+1} = CreatePlots(S, varargin{iplot}{1}, Sppt, varargin{iplot}{2:end});        
    end
end

if nargout >= 1,
    sOut = struct(...
        'listS', S ...
        ,'listHfig', {listHfig} ... % how to put a cell array in a struct field
        ,'Sppt', Sppt ...
        ,'probeh', probeh ...
        ,'rmsdDMv', rmsdDMv ...
        ,'itnum_min', itnum_min ...
        ,'fPlotNormIntensity', @PlotNormIntensity ...
        ,'fPlotBeta', @PlotBeta ...
        ,'fPlotProbeh', @PlotProbeh ...
        ,'fPlotRMSdDMv', @PlotRMSdDMv ...
        );
end

more on

end % main

function [hfig, hax, sCmetrics] = CreatePlots(S, sDisplayFun, Sppt, varargin)
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
                sCmetrics(ii) = sMtmp;
                
                if ~isempty(hfig)
                    figscale = CalcFigscale(hfig, figheight);
                    set(hfig, 'Position', figscale*get(hfig,'position'));
                    if ispc,
                        htmp = Sppt.CopyFigNewSlide(hfig);
                        %set(htmp,'Height',figheight);
                    else
                        saveas(hfig, [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
                    end
                end % if hfig
                
            end % for ii iter
        
            if strcmp(sMtmp(1).type, 'dEfields'),
                nw = S(1).NofW; % for convenience
                hfig_de = figure;
                hh = semilogy([S(2:N).iter], [sCmetrics.rmsdE_t].^2, '-', ...
                    [S(2:N).iter], mean([sCmetrics.rmsdE_t].^2, 1), '-', ...
                    [S(2:N).iter], [sCmetrics.rmsdE_m].^2, '--', ...
                    [S(2:N).iter], mean([sCmetrics.rmsdE_m].^2, 1), '--');
                
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
                
                newslide = Sppt.NewSlide(2);
                Sppt.CopyFigSlide(newslide, hfig_de);

            end
            
        case 'DisplayCEfields'

            for ii = 1:N-1,
                
                [hfig, hax, sCtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:},'hfig',hfig);
                if ~isempty(sCtmp),
                    sCmetrics(ii) = sCtmp;
                end
                
                if ~isempty(hfig),
                    figscale = CalcFigscale(hfig, figheight);
                    set(hfig, 'Position', figscale*get(hfig,'position'));
                    if ispc,
                        htmp = Sppt.CopyFigNewSlide(hfig);
                        %set(htmp,'Height',figheight);
                    else
                        saveas(hfig, [save_pn 'it' num2str(S(ii).iter) '_' sDisplayFun '.jpg']);
                    end % if ispc
                    
                end % if hfig
            end % for ii iter

            hfig_ce = figure;
            plotampphase([S(2:N).iter], [sCmetrics.CC],...
                'xlabel','Iteration #','title',[trialname ', \DeltaE Testbed Model Correlation (CC)']);
            newslide = Sppt.NewSlide(2);
            Sppt.CopyFigSlide(newslide, hfig_ce);
            
        otherwise, % one call per iteration
            sCmetrics = struct;
            for ii = 1:N,
                hfig = S(ii).(sDisplayFun)(varargin{:},'hfig',hfig);
                
                % check for valid figure
                if isempty(hfig),
                    continue
                end
                    
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

function [probeh, hfig, hax] = PlotProbeh(S, varargin)
     % [probeh, hfig, hax] = PlotProbeh(S, varargin)
     
     hfig = CheckOption('hfig', figure, varargin{:});
     hax = CheckOption('hax', [], varargin{:});

     % get texp
     [hfig, haxtexp, itnum_texp, texp] = PlotTexp(S);

     [itnum, probeh] = deal(zeros(size(S)));
     for ii = 1:length(S)
         itnum(ii)  = S(ii).iter;
         if ~isempty(S(ii).ImKeys),
             probeh(ii) = FitsGetKeywordVal(S(ii).ImKeys,'PROBEH');
         else
             probeh(ii) = NaN;
         end
     end

     figure(hfig);
     if ~isempty(hax), axes(hax); else, hax = gca; end
     yyaxis left
     semilogy(itnum, probeh, '-o'), grid on
     xlabel('Iteration #')
     ylabel('probeh')
     
     yyaxis right
     semilogy(itnum_texp, texp, '-x'), grid on
     ylabel('T_{exp} (s)')
     
end % PlotProbeh

function [hfig, hax, itnum, texp] = PlotTexp(S)

    [itnum, texp] = deal(zeros(size(S)));
    for ii = 1:length(S)
        itnum(ii) = S(ii).iter;
        if ~isempty(S(ii).ImKeys)
            texp(ii) = FitsGetKeywordVal(S(ii).ImKeys, 'texp');
        else
            % empty instance
            texp(ii) = NaN;
        end
    end
    
    hfig = figure;
    plot(itnum, texp, 'o'), grid
    xlabel('Iteration #')
    ylabel('T_{exp} (s)')
    hax = gca;
    
end % PlotTexp

function [betaused, betamin, hfig, hax] = PlotBeta(listS, varargin)
    % [betaused, betamin, hfig, hax] = PlotBeta(listS, varargin)
    
    hfig = CheckOption('hfig', figure, varargin{:});
    hax = CheckOption('hax', [], varargin{:});
    itnum = CheckOption('itnum', [listS.iter], varargin{:}); % use [listS.iter] - listS(1).iter to start with 0
    
    [ireg0, betaused, betamin] = deal(zeros(length(itnum),1));
    for ii = 1:length(listS)
        if ~isempty(listS(ii).ReducedKeys),
            ireg0(ii) = FitsGetKeywordVal(listS(ii).ReducedKeys, 'IREG0');
            if ireg0(ii) == -1,
                % optimal beta
                betaused(ii) = FitsGetKeywordVal(listS(ii).ReducedKeys, 'BMIN');
            else
                betaused(ii) = FitsGetKeywordVal(listS(ii).ReducedKeys, ['BSCAN' num2str(ireg0(ii), '%03d')]);
            end
            
            betamin(ii) = FitsGetKeywordVal(listS(ii).ReducedKeys, 'BMIN');
            
        else
            % empty instance
            betaused(ii) = NaN;
            betamin(ii) = NaN;
        end
    end % 
    
    figure(hfig);
    if ~isempty(hax), axes(hax); else, hax = gca; end
    hll = plot(itnum, betaused, '-or', itnum, betamin, '-xb');
    set(hll, 'LineWidth', 2);
    grid on
    xlabel('Iteration #')
    ylabel('Regularization \beta')
    legend('\beta used','\beta optimal')
    
end % PlotBeta

function [hfig, hax, han, itnum_min] = PlotNormIntensity(listS, varargin)

    hfig = CheckOption('hfig', figure, varargin{:});
    hax = CheckOption('hax', [], varargin{:});
    itnum = CheckOption('itnum', [listS.iter], varargin{:}); % use [listS.iter] - listS(1).iter to start with 0
    
    [NInt_co, NInt_inco, NInt_total] = deal(zeros(length(itnum), max([listS.Nlamcorr]) ));
    NInt_mean = zeros(length(itnum),1);
    for ii = 1:length(itnum)
        
        %         if isempty(listS(ii).NormIntensity_total)
        %             listS(ii).ReadImageCube;
        %         end
        %         if isempty(listS(ii).NormIntensity_co)
        %             listS(ii).ReadReducedCube;
        %         end
        %
        %
        %         for ilam = 1:listS(ii).Nlamcorr,
        %             NInt_co(ii, ilam) = listS(ii).NormIntensity_co(ilam);
        %             NInt_inco(ii, ilam) = listS(ii).NormIntensity_inco(ilam);
        %             NInt_total(ii, ilam) = listS(ii).NormIntensity_total(ilam);
        %         end % for ilam
        
        sC = listS(ii).GetContrast('display',false);
        if ~isempty(sC.co_lam_NI),
            NInt_co(ii, :) = sC.co_lam_NI;
            NInt_inco(ii, :) = sC.inco_lam_NI;
            NInt_total(ii, :) = sC.score_lam;
            NInt_mean(ii) = sC.mean;
        else
            NInt_co(ii, :) = NaN; % so it's not plotted
            NInt_inco(ii, :) = NaN;
            NInt_total(ii, :) = NaN;            
            NInt_mean(ii) = NaN;
        end
        
    end % ii

    figure(hfig);
    if ~isempty(hax), axes(hax); else, hax = gca; end
    hl = semilogy(itnum, NInt_total, '-', itnum, NInt_inco, '--', itnum, NInt_co, ':'); grid on
    xlabel('Iteration #')
    ylabel('Normalized Intensity')
    set(hl,'linewidth', 2)
    
    for ilam = 1:length(listS(1).lambda)
        legstr_total{ilam} = ['Total ' num2str(listS(1).lambda(ilam)/listS(1).NM, '%.0f') 'nm'];
        legstr_unmod{ilam} = ['Unmodulated ' num2str(listS(1).lambda(ilam)/listS(1).NM, '%.0f') 'nm'];
        legstr_mod{ilam}   = ['Modulated ' num2str(listS(1).lambda(ilam)/listS(1).NM, '%.0f') 'nm'];
    end
    %legend('Total', 'Unmodulated', 'Modulated')
    legend(legstr_total{:},legstr_unmod{:},legstr_mod{:})
    [itnum_min, NInt_total_min, NInt_inco_min, NInt_co_min, NInt_mean_min] = mindata(NInt_mean, itnum, mean(NInt_total, 2), mean(NInt_inco, 2), mean(NInt_co, 2), NInt_mean);
    han = FigureTitle(['Iter #' num2str(itnum_min) '; NI = ' num2str(NInt_mean_min,'%.1e') '; Mod = ' num2str(NInt_co_min,'%.1e ') '; Unmod = ' num2str(NInt_inco_min,'%.1e ')],'FontSize',12);

    
end % PlotNormIntensity

function [rmsdDMv, hfig, hax] = PlotRMSdDMv(listS, varargin)
    % [hfig, hax, rmsdDMv] = PlotRMSdDMv(listS, varargin)

    hfig = CheckOption('hfig', figure, varargin{:});
    hax = CheckOption('hax', [], varargin{:});
    itnum = CheckOption('itnum', [listS.iter], varargin{:}); % use [listS.iter] - listS(1).iter to start with 0

    % plotting differential, start with itnum(2)
    itnum_plot = itnum(2:end);
    itnum_plot = itnum_plot(:);
    
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
            if ~isempty(listS(ii).DMvCube),
                DMvtmp{idm} = squeeze(listS(ii).DMvCube{idm}(:,:,1));
            else
                % empty instance
                DMvtmp{idm} = NaN;
            end
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
    
    figure(hfig);
    if ~isempty(hax), axes(hax); else, hax = gca; end
    hh = semilogy(itnum_plot, rmsdDMv, '-o', itnum_plot, mean(rmsdDMv,2), '--');
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