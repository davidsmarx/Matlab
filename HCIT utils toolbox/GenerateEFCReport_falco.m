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
Sppt = CheckOption('Sppt', 'new', varargin{:});
bSpptIsNew = ischar(Sppt) && strcmpi(Sppt, 'new'); % Sppt gets overwritten below, save this before that happens
run_bn = CheckOption('run_bn', ['falco_testbed_run' num2str(runnum)], varargin{:});
listSin = CheckOption('listS', [], varargin{:}); % if listS of CfalcoRunData for iterations already exists

% % check if PowerPoint Presentation already exists and still there
% try
%     Sppt.Presentation.Slides,
% catch
%     clear Sppt; Sppt = [];
% end

% open PowerPoint if necessary, and plots are requested (Windows pc only)
if bSpptIsNew && ispc %&& ~isempty(varargin),
    Sppt = Cppt(ppt_fn);
end

% initial mp is empty, gets read by first iteration
mp = [];

% load all the falco data for the iterations, if not input
if isempty(listSin)
    S = LoadIterations(runnum, TrialNum, listItnum, mp, varargin{:});
else
    % check that S contains the requested iterations
    [listItnum, i_S_list_iter] = intersect([listSin.iter], listItnum);
    if isempty(listItnum)
        error('requested listItnum not in provided iterations');
    end
    S = listSin(i_S_list_iter);
end

% some path definitions
report_pn = PathTranslator(fullfile(getenv("DATA_ROOT"), run_bn, 'reports', S(1).runLabel));

% plot graphs of metrics v itnum on the first slide
if ~isempty(Sppt), slide = Sppt.NewSlide(1); end
hfig = figure_mxn(2,2);
set(hfig, 'position', [322 263 1800 1000]);
hax(1,1) = subplot(2,2,1); [~, ~, ~, itnum_min] = PlotNormIntensity(S, 'hfig', hfig, 'hax', hax(1,1));
hax(1,2) = subplot(2,2,2); PlotBeta(S, 'hfig', hfig, 'hax', hax(1,2));
hax(2,1) = subplot(2,2,3); probeh = PlotProbeh(S, 'hfig', hfig, 'hax', hax(2,1));
hax(2,2) = subplot(2,2,4); rmsdDMv = PlotRMSdDMv(S, 'hfig', hfig, 'hax', hax(2,2));
if ~isempty(Sppt),
    hPic = Sppt.CopyFigSlide(slide, hfig);
else
    fSaveas(hfig, report_pn, 'summary', ['summary_it' num2str(S(1).iter) '_it' num2str(S(end).iter)], []);
end

% compare full-band (mean over subbands/modes) normalized intensity images
% -- Total (unprobed), Modulated, Unmodulated -- side by side between two
% iterations (default: first and last iteration of this report)
itnum_compare1 = CheckOption('itnum_compare1', S(end).iter-8, varargin{:});
itnum_compare2 = CheckOption('itnum_compare2', S(end).iter-4, varargin{:});
[hfig_cmp, ~, sImageCompareData] = PlotCompareNormInt(S, itnum_compare1, itnum_compare2);
if ~isempty(Sppt)
    newslide = Sppt.NewSlide([]); % append at end
    Sppt.CopyFigSlide(newslide, hfig_cmp);
end
if ~exist(report_pn, 'dir'), mkdir(report_pn); end
save(fullfile(report_pn, ['compare_normint_it' num2str(itnum_compare1) '_it' num2str(itnum_compare2) '.mat']), '-struct', 'sImageCompareData');

% add saved falco figures
% only when starting a new PowerPoint (option 'Sppt' == 'new');
% skip if continuing/appending to an already-open presentation passed in by the caller
figures_pn = [S(1).Rundir_pn '/figures'];
if bSpptIsNew && ~isempty(Sppt) && exist(PathTranslator(figures_pn), 'dir')
    listPng = dir(fullfile(PathTranslator(figures_pn), '*.png'));
    for ii = 1:length(listPng)
        slide = Sppt.NewSlide(1+ii);
        fn = fullfile(listPng(ii).folder, listPng(ii).name);
        hh = invoke(slide.Shapes, 'AddPicture', fn, true, true, 100, 100);
        % hh.Left, hh.Top, hh.Width
    end
end

% call the plotting methods
listHfig = {};
for iplot = 1:length(varargin),
    if iscell(varargin{iplot}),
        listHfig{end+1} = CreatePlots(S, varargin{iplot}{1}, Sppt, varargin{iplot}{2:end}, 'save_pn', report_pn);        
    end
end

% save Sppt
if ~isempty(Sppt)
    fn = PathTranslator(fullfile(getenv("DATA_ROOT"), run_bn, 'reports', [S(1).runLabel '_it' num2str(S(1).iter) '_' num2str(S(end).iter) '.pptx']));

    % Check if file already exists
    if exist(fn, 'file')
        % Prompt user for action
        choice = questdlg(['File already exists: ' fn], ...
            'File Exists', ...
            'Replace', 'Save with new name', 'Cancel', 'Cancel');

        switch choice
            case 'Replace'
                % Delete existing file and proceed with save
                delete(fn);
                Sppt.Presentation.SaveAs(fn);
                fprintf('Replaced existing file: %s\n', fn);

            case 'Save with new name'
                % Generate new filename with datetime
                [fpath, fname, fext] = fileparts(fn);
                datetime_str = datestr(now, 'yyyymmdd_HHMMSS');
                fn_new = fullfile(fpath, [fname '_' datetime_str fext]);
                Sppt.Presentation.SaveAs(fn_new);
                fprintf('Saved as: %s\n', fn_new);

            case 'Cancel'
                % Do nothing
                fprintf('Save cancelled. File not saved.\n');

            otherwise
                % User closed dialog - do nothing
                fprintf('Save cancelled. File not saved.\n');
        end
    else
        % File doesn't exist, save normally
        Sppt.Presentation.SaveAs(fn);
        fprintf('Saved: %s\n', fn);
    end
end

if nargout >= 1,
    sOut = struct(...
        'listS', S ...
        ,'listHfig', {listHfig} ... % how to put a cell array in a struct field
        ,'Sppt', Sppt ...
        ,'probeh', {probeh} ...
        ,'rmsdDMv', {rmsdDMv} ...
        ,'itnum_min', itnum_min ...
        ,'fPlotNormIntensity', @PlotNormIntensity ...
        ,'fPlotBeta', @PlotBeta ...
        ,'fPlotProbeh', @PlotProbeh ...
        ,'fPlotRMSdDMv', @PlotRMSdDMv ...
        ,'fPlotCompareNormInt', @PlotCompareNormInt ...
        );
end

more on

end % main

function S = LoadIterations(runnum, TrialNum, listItnum, mp, varargin)
%

% options
max_empties = CheckOption('max_empties', 3, varargin{:});

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
        S(ii) = CfalcoRunData(runnum, TrialNum, listItnum(ii), 'mp', mp, varargin{:});

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

    % if last iteration is empty, remove it
    if isempty(S(end).ImCube)
        S = S(1:end-1);
    end
    
elseif isa(listItnum, 'CfalcoRunData')
    S = listItnum;
else
    error(['listItnum type error: ' class(listItnum)]);
end

end % function LoadIterations

function [hfig, hax, sCmetrics] = CreatePlots(S, sDisplayFun, Sppt, varargin)
    % create the plots
    % some plots are differential
    % some plots we also plot metrics v itnum

    save_pn = CheckOption('save_pn', pwd, varargin{:}); % if ~ispc, must be absolute path
    figheight = CheckOption('figheight', 700, varargin{:}); % for ppt display
    trialname = CheckOption('trialname', '', varargin{:});

    % create path to put plots, if necessary
    if ~ispc && ~exist(save_pn)
        mkdir(save_pn)
        if isunix
            % set permission
            command = ['chmod 775 ' save_pn];
            system(command)
        end
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
        case 'DisplayDEfields'
            %if any(strcmp(sDisplayFun, listDiff)),
            for ii = 1:N-1,
                % each iteration, create a new DisplayDEfields
                
                [hfig, hax, sMtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:},'hfig',hfig);
                if ~isempty(sMtmp)
                    sCmetrics(ii) = sMtmp;
                end
                
                if ~isempty(hfig)
                    figscale = CalcFigscale(hfig, figheight);
                    set(hfig, 'Position', figscale*get(hfig,'position'));

                    % if using ImageCube to toggle model/measure, make the gif
                    if ~isempty(get(hfig, 'KeyPressFcn'))
                        % send the 'G' keystroke to create the gif
                        fungif = get(hfig, 'KeyPressFcn');
                        
                        % save_pn must be full path
                        gif_fn = fullfile(save_pn, sDisplayFun, ['it_' num2str(S(ii).iter) '.gif']);
                        pn = fileparts(gif_fn);
                        if ~exist(pn, 'dir'),
                            mkdir(pn);
                            if isunix
                                % set permission
                                command = ['chmod 775 ' pn];
                                system(command)
                            end

                        end
                        
                        % this sends an event to hfig to create the gif
                        thisevent = struct('Modifier', 'shift', 'Key', 'g', 'gif_fn', gif_fn);
                        fungif(hfig, thisevent);
                        
                        % check and insert gif to PowerPoint
                        if exist(gif_fn, "file") && ispc
                            hh = Sppt.AddPictureNewSlide(gif_fn);
                        end % if copy to PowerPoint

                    else
                        fSaveas(hfig, save_pn, sDisplayFun, 'it', S(ii).iter);

                    end % if ImageCube

                end % if hfig
                
            end % for ii iter
        
            % plot a summary rms dE vs iteration
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

        case 'DisplayDMv'
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
        
        case 'DisplayCEfields'

            for ii = 1:N-1,
                try
                    [hfig, hax, sCtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:},'hfig',hfig);
                catch ME
                    hfig = [];                    
                    disp(ME.message);
                    for ii = 1:length(ME.stack),
                        disp(ME.stack(ii));
                    end
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
    % fSaveas(hfig, save_pn, sDisplayFun, bn, iter)
    %
    % save figure to specified location as .jpg and .fig

    fn = fullfile(save_pn, sDisplayFun, [bn '_' num2str(iter) '.jpg']);
    fnfig = fullfile(save_pn, sDisplayFun, [bn '_' num2str(iter) '.fig']);

    pn = fileparts(fn);
    if ~exist(pn, 'dir'),
        mkdir(pn);
        if isunix
            % set permission
            command = ['chmod 775 ' pn];
            system(command)
        end
    end
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
     probeh = cell(1, S(1).NofW); % zeros(length(S),
     for ii = 1:length(S)
         for iw = 1:S(ii).NofW
             for ip = 1:S(ii).Nppair
                 Itmp(:,ip) = S(ii).ProbeMeasAmp{iw, ip}(S(ii).bMask).^2;
             end % each probe
             probeh{iw}(ii) = mean(Itmp(:));
         end % each subband
         itnum(ii) = S(ii).iter;
     end % each iteration

     figure(hfig);
     if ~isempty(hax), axes(hax); else, hax = gca; end
     yyaxis left
     % S.Nlamcorr = mp.Nsbp;
     % S.NofW = S.Nlamcorr * S.Nstar
     % imode = (istar-1)*S.Nlamcorr + iwv;
     strStar = {'On-axis', 'Off-axis'};
     symbols = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', '+', 'x'};     
     for imode = 1:length(probeh)
         ibnd = rem(imode-1, S(1).Nlamcorr) + 1;
         istar = floor((imode - ibnd)./S(1).Nlamcorr) + 1;
         hl(imode) = semilogy(itnum, probeh{imode}, ['-' symbols{imode}]);
         grid on, hold on
         legstr{imode} = ['Band ' num2str(ibnd) ', ' strStar{istar} ' Star'];
     end
     %set(gca, 'YScale', 'log')
     hold off
     legend(legstr{:})
     xlabel('Iteration #')
     ylabel('Mean Probe Intensity')
     
     yyaxis right
     semilogy(itnum_texp, texp, '-x'), grid on
     ylabel('T_{exp} (s)')
     legstr{end+1} = 'T_{exp}';

     legend(legstr{:})
      
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

    % look for falcoData.ctrl.log10regHist and falcoData.ctrl.dmfacHist
    betaused = listS(end).falcoData.log10regHist(itnum);
    betamin = [];
    dmfac_used = listS(end).falcoData.ctrl.dmfacHist(itnum);
    
    if isempty(hfig),
        hfig = figure;
    else, 
        figure(hfig);
    end
    if ~isempty(hax), axes(hax); else, hax = gca; end

    % left axis is beta
    yyaxis left
    hll = plot(itnum, betaused, '-o');
    set(hll, 'LineWidth', 2);
    grid on
    xlabel('Iteration #')
    ylabel('Regularization \beta')

    % right axis is ddm gain
    yyaxis right
    hrr = plot(itnum, dmfac_used, '-x');
    set(gca,'YScale','log')
    set(hrr, 'LineWidth', 2);
    ylabel('\deltaDM Scale Factor')

    legend('\beta', '\deltaDM Scale')

    
end % PlotBeta

function [hfig, hax, han, itnum_min] = PlotNormIntensity(listS, varargin)
    % [hfig, hax, han, itnum_min] = PlotNormIntensity(listS, varargin)
    % CheckOption('hfig', [], varargin{:});
    % CheckOption('hax', [], varargin{:});
    % CheckOption('itnum', [listS.iter], varargin{:}); % use [listS.iter] - listS(1).iter to start with 0
    % CheckOption('ylim', [], varargin{:});

    hfig = CheckOption('hfig', [], varargin{:});
    hax = CheckOption('hax', [], varargin{:});
    itnum = CheckOption('itnum', [listS.iter], varargin{:}); % use [listS.iter] - listS(1).iter to start with 0
    ylim = CheckOption('ylim', [], varargin{:});

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
    if isfield(listS(1).mp, 'toggledMSWC') && listS(1).mp.toggledMSWC
        NInt_mean = sum(NInt_total, 2); % mean of subbands and stars. Not correct for toggled
    else
        NInt_mean = mean(NInt_total, 2); % mean of subbands and stars. Not correct for toggled
    end
    
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

    % auto ylim
    if isempty(ylim)
        ylim = [1e-10 0] + [0 1].*get(hax, 'ylim');
    end
    set(hax, 'ylim', ylim)

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

function [hfig, hax, sImageData] = PlotCompareNormInt(S, itnum1, itnum2, varargin)
    % [hfig, hax, sImageData] = PlotCompareNormInt(S, itnum1, itnum2, varargin)
    %
    % compare full-band (mean over subbands/modes) normalized intensity
    % images between two iterations: 3 rows x 2 cols
    %    row 1 = Total (unprobed)
    %    row 2 = Modulated
    %    row 3 = Unmodulated
    %    col 1 = itnum1, col 2 = itnum2
    %
    % S = array of CfalcoRunData objects (already loaded)
    % sImageData = struct with the images used in the plot, for saving to .mat
    %
    % Options:
    %
    %     CheckOption('bLog', true, varargin{:});
    %     CheckOption('clim', [-9 -6.5], varargin{:});
    %     CheckOption('xlim', S(1).XlimDefault, varargin{:});
    %     CheckOption('ylim', S(1).YlimDefault, varargin{:});
    %     CheckOption('title_left', ['Iter #' num2str(itnum1)])
    %     CheckOption('title_right', ['Iter #' num2str(itnum2)])
    %     CheckOption('save_fn', [], varargin{:});

    bLog = CheckOption('bLog', true, varargin{:});
    climopt = CheckOption('clim', [-9 -6.5], varargin{:});
    xlim = CheckOption('xlim', S(1).XlimDefault, varargin{:});
    ylim = CheckOption('ylim', S(1).YlimDefault, varargin{:});
    title1 = CheckOption('title_left', ['Iter #' num2str(itnum1)], varargin{:});
    title2 = CheckOption('title_right', ['Iter #' num2str(itnum2)], varargin{:});
    save_fn = CheckOption('save_fn', [], varargin{:});

    i1 = find([S.iter] == itnum1, 1);
    i2 = find([S.iter] == itnum2, 1);
    if isempty(i1) || isempty(i2)
        error('PlotCompareNormInt: itnum1 (%d) or itnum2 (%d) not found in S', itnum1, itnum2);
    end
    Scol = [S(i1) S(i2)];
    itnumcol = [itnum1 itnum2];
    Ncols = length(itnumcol);
    titlecol = {title1, title2};

    [x, y] = CreateGrid(Scol(1).ImCubeUnProbFullBand, 1./Scol(1).ppl0);

    rowprop  = {'ImCubeUnProbFullBand', 'CohIntFullBand', 'IncIntFullBand'};
    rowlabel = {'Total (UnProbed)', 'Modulated', 'Unmodulated'};
    rowtitle = {'Total Norm Intensity', 'Modulated', 'Unmodulated'};

    hfig = figure_mxn(3, Ncols);
    hax = zeros(3, Ncols);
    %sImageData = struct('itnum1', itnum1, 'itnum2', itnum2, 'x', x, 'y', y);
    sImageData = struct;

    % explicit layout so rows are packed tightly, with extra room at the
    % top of each column reserved for a bold supertitle (e.g. 'On-Axis Star')
    left_margin   = 0.08;
    right_margin  = 0.13;
    hgap          = 0.05;
    top_margin    = 0.075;
    bottom_margin = 0.09;
    vgap          = 0.05;

    col_width  = (1 - left_margin - right_margin - (Ncols-1)*hgap) / Ncols;
    row_height = (1 - top_margin - bottom_margin - 2*vgap) / 3;

    for icol = 1:Ncols
        sImageData(icol).itnum = itnumcol(icol);
        sImageData(icol).x = x;
        sImageData(icol).y = y;

        left = left_margin + (icol-1)*(col_width + hgap);

        for irow = 1:3
            bottom = 1 - top_margin - irow*row_height - (irow-1)*vgap;
            hax(irow,icol) = subplot('Position', [left bottom col_width row_height]);
            im = Scol(icol).(rowprop{irow});
            sImageData(icol).(rowprop{irow}) = im;

            if bLog
                imageschcit(x, y, log10(abs(im))); axis image
            else
                imageschcit(x, y, im); axis image
            end
            if ~isempty(climopt), set(gca,'clim',climopt); end
            if ~isempty(xlim), set(gca,'xlim',xlim); end
            if ~isempty(ylim), set(gca,'ylim',ylim); end

            if icol == 1
                ylabel('\lambda/D')
            end
            if irow == 3
                xlabel('\lambda/D')
            end

            title(rowtitle{irow})

            % only the last column shows a colorbar; add it without
            % shrinking that column's axes so all columns stay the same size
            axpos = get(gca, 'Position');
            if icol == Ncols
                hcb = colorbar;
                colorbartitle(hcb, 'log_{10} Norm Intensity');
                set(gca, 'Position', axpos);
                set(hcb, 'Position', [axpos(1)+axpos(3)+0.015, axpos(2), 0.02, axpos(4)]);
            end

        end % for each column (iteration)

        % bold, obvious supertitle for this column, e.g. 'On-Axis Star'
        annotation(hfig, 'textbox', [left, 1-0.06, col_width, 0.05] ...
            , 'String', titlecol{icol} ...
            , 'FontSize', 24 ...
            , 'FontWeight', 'bold' ...
            , 'Color', 'r' ...
            , 'HorizontalAlignment', 'center' ...
            , 'VerticalAlignment', 'middle' ...
            , 'EdgeColor', 'none' ...
            );

        % row label
        ylpos = get(get(hax(irow,1),'YLabel'),'Position');
        text(hax(irow,1), ylpos(1) - 2, ylpos(2), rowlabel{irow} ...
            , 'Rotation', 90 ...
            , 'HorizontalAlignment', 'center' ...
            , 'VerticalAlignment', 'bottom' ...
            , 'FontSize', 14 ...
            , 'Color', 'b' ...
            , 'FontWeight', 'bold' ...
            );
    end % for each row (norm intensity type)

    if ~isempty(save_fn)
        fprintf('Saving image data to %s...', save_fn);
        save(save_fn, "sImageData");
    end
    
end % PlotCompareNormInt