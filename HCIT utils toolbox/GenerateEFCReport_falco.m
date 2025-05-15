function sOut = GenerateEFCReport_falco(runnum, TrialNum, listItnum, varargin)
% sOut = GenerateEFCReport(runnum, listItnum, cDisplayFun)
%
% listItnum can be array of itnums, or array of CfalcoRunData objects
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
Sppt = CheckOption('Sppt', [], varargin{:});
max_empties = CheckOption('max_empties', 3, varargin{:});
run_pn = CheckOption('run_pn', 'falco_testbed_run', varargin{:});

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

% initial mp is empty, gets read by first iteration
mp = [];

% get the CfalcoRunData objects
% if 3 successive iterations have no data, stop
cnt_empty = 0;
N = length(listItnum);
if isnumeric(listItnum),
    %if isscalar(listItnum)
    %ii = 1;
    %while true,
    for ii = 1:N;
        fprintf('reading itnum %d\n', listItnum(ii));
        S(ii) = CfalcoRunData(runnum, TrialNum, listItnum(ii), 'mp', mp, 'run_pn', run_pn);

        if isempty(S(ii).ImCube)
            cnt_empty = cnt_empty + 1;
        else
            % reset, only count successive empties
            cnt_empty = 0;
        end
        
        if cnt_empty >= max_empties
            S = S(1:end-max_empties);
            break
        end
        
        % update mp, so it is used for next iteration and avoid reading
        % config every iteration
        mp = S(ii).mp;
        
        % check if this is the last iteration
        % out.Itr from "_snippet.mat" is last iteration
        if listItnum(ii) >= S(ii).falcoData.Itr, 
            break;
        end
        
    end
elseif isa(listItnum, 'CfalcoRunData')
    S = listItnum;
else
    error(['listItnum type error: ' class(listItnum)]);
end

%
%saveas_pn = ['./falco_testbed_run' num2str(S(1).runnum) '/data/' S(1).runLabel '/figures'];
saveas_pn = ['./' run_pn num2str(S(1).runnum) '/data/' S(1).runLabel '/figures'];

% plot graphs of metrics v itnum on the first slide
if ~isempty(Sppt), slide = Sppt.NewSlide(1); end
hfig = figure_mxn(2,2);
hax(1,1) = subplot(2,2,1); [~, ~, ~, itnum_min] = PlotNormIntensity(S, 'hfig', hfig, 'hax', hax(1,1));
hax(1,2) = subplot(2,2,2); PlotBeta(S, 'hfig', hfig, 'hax', hax(1,2));
hax(2,1) = subplot(2,2,3); probeh = PlotProbeh(S, 'hfig', hfig, 'hax', hax(2,1));
hax(2,2) = subplot(2,2,4); rmsdDMv = PlotRMSdDMv(S, 'hfig', hfig, 'hax', hax(2,2));
if ~isempty(Sppt),
    hPic = Sppt.CopyFigSlide(slide, hfig);
else
    fSaveas(hfig, saveas_pn, 'summary', ['summary_it' num2str(S(1).iter) '_it' num2str(S(end).iter)], []);
end

% add saved falco figures
list_fignum_to_copy = [1 2 51 91 401];
figures_pn = [S(1).Rundir_pn '/figures'];
if exist(PathTranslator(figures_pn), 'dir')
    %listPng = dir(PathTranslator([figures_pn '/*.png']));
    for ii = 1:length(list_fignum_to_copy) %length(listPng)
        if ~isempty(Sppt), slide = Sppt.NewSlide(1+ii); end
        fn = fullfile(PathTranslator(figures_pn), ['figure_' num2str(list_fignum_to_copy(ii)) '.png']);
        if ~exist(fn, 'file')
            continue
        end
        if ~isempty(Sppt)
            hh = invoke(slide.Shapes, 'AddPicture', fn, true, true, 100, 100);
            % hh.Left, hh.Top, hh.Width
        end
    end
end

% call the plotting methods
listHfig = {};
for iplot = 1:length(varargin),
    if iscell(varargin{iplot}),
        listHfig{end+1} = CreatePlots(S, varargin{iplot}{1}, Sppt, varargin{iplot}{2:end}, 'save_pn', saveas_pn);        
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

    save_pn = CheckOption('save_pn', ['./' sDisplayFun '/'], varargin{:}); % if ~ispc
    figheight = CheckOption('figheight', 700, varargin{:}); % for ppt display
    trialname = CheckOption('trialname', '', varargin{:});

    % create path to put plots, if necessary
    if ~ispc && ~exist(save_pn)
        mkdir(save_pn)
    end

    % list of CfalcoRunData methods where the first argument is a reference
    % iteration:
    listDiff = {
        'DisplayDEfields'
        'DisplayDMv'
        };
    
    N = length(S);

    % if only 1 iteration, can't do displays that use differences
    if N <= 1 && any(strcmp(sDisplayFun, [listDiff(:); {'DisplayCEfields'}]))
        % just return
        hfig = []; hax = []; sCmetrics = struct;
        return
    end
    
    hfig = [];
    switch sDisplayFun,
        case listDiff
            %if any(strcmp(sDisplayFun, listDiff)),
            for ii = 1:N-1,
                [hfig, hax, sMtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:},'hfig',hfig);
                if ~isempty(sMtmp)
                    sCmetrics(ii) = sMtmp;
                end
                
                if ~isempty(hfig)
                    figscale = CalcFigscale(hfig, figheight);
                    set(hfig, 'Position', figscale*get(hfig,'position'));
                    if ispc,
                        htmp = Sppt.CopyFigNewSlide(hfig);
                        %set(htmp,'Height',figheight);
                    else
                        fSaveas(hfig, save_pn, sDisplayFun, 'it', S(ii).iter);
                    end
                end % if hfig
                
            end % for ii iter
        
            if any(strcmp({sCmetrics.type}, 'dEfields')),
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
                
                if ispc
                    newslide = Sppt.NewSlide(2);
                    Sppt.CopyFigSlide(newslide, hfig_de);
                else
                    fSaveas(hfig_de, save_pn, 'summary', ['magdE_it' num2str(S(1).iter) '_it' num2str(S(end).iter)], []);
                end
            end
            
        case 'DisplayCEfields'

            for ii = 1:N-1,
                try
                    [hfig, hax, sCtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:},'hfig',hfig);
                catch ME
                    hfig = [];                    
                    disp(ME.message);
                    disp(ME.stack);
                end

                if ~isempty(sCtmp),
                    sCmetrics(ii) = sCtmp;
                end
                
                if isempty(hfig)
                    continue
                end

                figscale = CalcFigscale(hfig, figheight);
                set(hfig, 'Position', figscale*get(hfig,'position'));
                if ispc,
                    htmp = Sppt.CopyFigNewSlide(hfig);
                    %set(htmp,'Height',figheight);
                else
                    fSaveas(hfig, save_pn, sDisplayFun, 'it', S(ii).iter);
                end % if ispc
                    
            end % for ii iter

            hfig_ce = figure;
            plotampphase([S(2:N).iter], [sCmetrics.CC],...
                'xlabel','Iteration #','title',[trialname ', \DeltaE Testbed Model Correlation (CC)']);
            
            if ~isempty(Sppt)
                newslide = Sppt.NewSlide(2);
                Sppt.CopyFigSlide(newslide, hfig_ce);
            else
                fSaveas(hfig_ce, save_pn, 'summary', ['CE_it' num2str(S(1).iter) '_it' num2str(S(end).iter)], []);
            end
            
        otherwise, % one call per iteration
            sCmetrics = struct;
            for ii = 1:N,
                try
                    hfig = S(ii).(sDisplayFun)(varargin{:},'hfig',hfig);
                catch ME
                    hfig = [];
                    disp(ME.message);
                    disp(ME.stack);
                end
                
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
                    fSaveas(hfig, save_pn, sDisplayFun, 'it', S(ii).iter);
                end
            end
            
    end % switch

end % CreatePlots

function fSaveas(hfig, save_pn, sDisplayFun, bn, iter)

    fn = fullfile(save_pn, sDisplayFun, [bn '_' num2str(iter) '.jpg']);
    fnfig = fullfile(save_pn, sDisplayFun, [bn '_' num2str(iter) '.fig']);

    pn = fileparts(fn);
    if ~exist(pn, 'dir'), mkdir(pn); end
    saveas(hfig, fn);
    saveas(hfig, fnfig);


end

function figscale = CalcFigscale(hfig, figheight)
    
    pos = get(hfig,'Position');
    ysize = pos(end);
    figscale = figheight/ysize;

end % CalcFigscale

function [probeh, hfig, hax] = PlotProbeh(S, varargin)
     % [probeh, hfig, hax] = PlotProbeh(S, varargin)
     
     hax = CheckOption('hax', [], varargin{:});

     % get texp
     if ~isempty(hax),
         [~, ~, itnum_texp, texp] = PlotTexp(S, 'nodisplay', true);
         hfig = hax.Parent;
     else
         [hfig, hax, itnum_texp, texp] = PlotTexp(S);
     end

     itnum = zeros(size(S));
     probeh = zeros(length(S), S(1).NofW);
     for ii = 1:length(S)
         for iw = 1:S(ii).NofW
             for ip = 1:S(ii).Nppair
                 Itmp(:,ip) = S(ii).ProbeMeasAmp{iw, ip}(S(ii).bMask).^2;
             end % each probe
             probeh(ii, iw) = mean(Itmp(:));
         end % each subband
         itnum(ii) = S(ii).iter;
     end % each iteration

     figure(hfig);
     if ~isempty(hax), axes(hax); else, hax = gca; end
     yyaxis left
     semilogy(itnum, probeh, '-o'), grid on
     xlabel('Iteration #')
     ylabel('Mean Probe Intensity')
     
     yyaxis right
     semilogy(itnum_texp, texp, '-x'), grid on
     ylabel('T_{exp} (s)')
     
end % PlotProbeh

function [hfig, hax, itnum, texp] = PlotTexp(S, varargin)

    nodisplay = CheckOption('nodisplay', false, varargin{:}); % in case you only want the texp data
    
    [itnum, texp] = deal(zeros(size(S)));
    for ii = 1:length(S)
        itnum(ii) = S(ii).iter;
        if ~isempty(S(ii).ReducedKeys)
            try
                texp(ii) = FitsGetKeywordVal(S(ii).ReducedKeys, 'texp1');
            catch
                texp(ii) = 0;
            end
        else
            % empty instance
            texp(ii) = NaN;
        end
    end
    
    if nodisplay
        hfig = []; hax = [];
        
    else
        hfig = figure;
        plot(itnum, texp, 'o'), grid
        xlabel('Iteration #')
        ylabel('T_{exp} (s)')
        hax = gca;
    end
    
end % PlotTexp

function [betaused, betamin, hfig, hax] = PlotBeta(listS, varargin)
    % [betaused, betamin, hfig, hax] = PlotBeta(listS, varargin)
    %
    % falco does not record betamin
    
    hfig = CheckOption('hfig', [], varargin{:});
    hax = CheckOption('hax', [], varargin{:});
    itnum = CheckOption('itnum', [listS.iter], varargin{:}); % use [listS.iter] - listS(1).iter to start with 0
    
    betaused = listS(end).falcoData.log10regHist(itnum);
    betamin = [];
    
    if isempty(hfig),
        hfig = figure;
    else, 
        figure(hfig);
    end
    if ~isempty(hax), axes(hax); else, hax = gca; end
    hll = plot(itnum, betaused, '-or');
    set(hll, 'LineWidth', 2);
    grid on
    xlabel('Iteration #')
    ylabel('Regularization \beta')
    legend('\beta used') %,'\beta optimal')
    
end % PlotBeta

function [hfig, hax, han, itnum_min] = PlotNormIntensity(listS, varargin)

    hfig = CheckOption('hfig', [], varargin{:});
    hax = CheckOption('hax', [], varargin{:});
    itnum = CheckOption('itnum', [listS.iter], varargin{:}); % use [listS.iter] - listS(1).iter to start with 0
    
    itnum = itnum(:); % force column vector
    [NInt_co, NInt_inco, NInt_total] = deal(zeros(length(itnum), max([listS.Nlamcorr]) ));
    NInt_mean = zeros(length(itnum),1);
    
    % %%% this parts calculates contrast or NI for each iteration
    % %%% not working right now, probably because of bMask problems
    % %%% use falcoData instead
    %     for ii = 1:length(itnum)
    %
    %         %         if isempty(listS(ii).NormIntensity_total)
    %         %             listS(ii).ReadImageCube;
    %         %         end
    %         %         if isempty(listS(ii).NormIntensity_co)
    %         %             listS(ii).ReadReducedCube;
    %         %         end
    %         %
    %         %
    %         %         for ilam = 1:listS(ii).Nlamcorr,
    %         %             NInt_co(ii, ilam) = listS(ii).NormIntensity_co(ilam);
    %         %             NInt_inco(ii, ilam) = listS(ii).NormIntensity_inco(ilam);
    %         %             NInt_total(ii, ilam) = listS(ii).NormIntensity_total(ilam);
    %         %         end % for ilam
    %
    %         %sC = listS(ii).GetContrast('display',false);
    %         if ~isempty(sC.co_lam_NI),
    %             NInt_co(ii, :) = sC.co_lam_NI;
    %             NInt_inco(ii, :) = sC.inco_lam_NI;
    %             NInt_total(ii, :) = sC.score_lam;
    %             NInt_mean(ii) = sC.mean;
    %         else
    %             NInt_co(ii, :) = NaN; % so it's not plotted
    %             NInt_inco(ii, :) = NaN;
    %             NInt_total(ii, :) = NaN;
    %             NInt_mean(ii) = NaN;
    %         end
    %
    %     end % ii

    NInt_co = listS(1).falcoData.normIntModScore(itnum, :); % column per band*star = listS(1).NofW
    NInt_inco = listS(1).falcoData.normIntUnmodScore(itnum, :); % column per band*star = listS(1).NofW
    NInt_total = listS(1).falcoData.normIntMeasScore(itnum, :); % column per band*star
    NInt_mean = mean(NInt_total, 2); % mean across the band ???
    
    if isempty(hfig),
        hfig = figure;
    else
        figure(hfig);
    end
    if ~isempty(hax), axes(hax); else, hax = gca; end
    hl = semilogy(itnum(:), NInt_total, '-', itnum, NInt_inco, '--', itnum, NInt_co, ':'); grid on
    xlabel('Iteration #')
    ylabel('Normalized Intensity')
    set(hl,'linewidth', 2)
    
    % S.Nlamcorr = mp.Nsbp;
    % S.NofW = S.Nlamcorr * S.Nstar
    strStar = {'On-axis', 'Off-axis'};
    for ibnd = 1:listS(1).NofW % check length of .lambda always equals Nsbd * Nstar ???
        % when there are multiple subbands with multiple stars, order will
        % be determined in falco routine. For now, one subband, two stars
        if listS(1).NofW == listS(1).Nstar,
            istar = ibnd;
            ilam = 1;
        else
            istar = 1;
            ilam = ibnd;
        end

        legstr_total{ibnd} = ['Total ' strStar{istar} ' ' num2str(listS(1).lambda(ilam)/listS(1).NM, '%.0f') 'nm'];
        legstr_unmod{ibnd} = ['Unmodulated ' strStar{istar} ' ' num2str(listS(1).lambda(ilam)/listS(1).NM, '%.0f') 'nm'];
        legstr_mod{ibnd}   = ['Modulated ' strStar{istar} ' ' num2str(listS(1).lambda(ilam)/listS(1).NM, '%.0f') 'nm'];
    end
    %legend('Total', 'Unmodulated', 'Modulated')
    legend(legstr_total{:},legstr_unmod{:},legstr_mod{:})
    [itnum_min, NInt_total_min, NInt_inco_min, NInt_co_min, NInt_mean_min] = mindata(NInt_mean, itnum, mean(NInt_total, 2), mean(NInt_inco, 2), mean(NInt_co, 2), NInt_mean);
    han = FigureTitle(['Trial # ' num2str(listS(1).trialNum) '; Iter #' num2str(itnum_min) '; NI = ' num2str(NInt_mean_min,'%.1e') '; Mod = ' num2str(NInt_co_min,'%.1e ') '; Unmod = ' num2str(NInt_inco_min,'%.1e ')],'FontSize',12);

    
end % PlotNormIntensity

function [rmsdDMv, hfig, hax] = PlotRMSdDMv(listS, varargin)
    % [hfig, hax, rmsdDMv] = PlotRMSdDMv(listS, varargin)

    hfig = CheckOption('hfig', [], varargin{:});
    hax = CheckOption('hax', [], varargin{:});
    itnum = CheckOption('itnum', [listS.iter], varargin{:}); % use [listS.iter] - listS(1).iter to start with 0

    % plotting differential, start with itnum(2)
    itnum_plot = itnum(2:end);
    itnum_plot = itnum_plot(:);
        
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
    
    if isempty(hfig),
        hfig = figure;
    else
        figure(hfig);
    end
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