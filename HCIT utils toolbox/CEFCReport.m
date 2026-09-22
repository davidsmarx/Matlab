classdef CEFCReport < handle
    % CEFCReport - EFC iteration report generator
    %
    % R = CEFCReport(runnum, TrialNum, varargin)
    %
    % Handle-class replacement for the free-function GenerateEFCReport_falco.m.
    % The constructor auto-discovers every valid iteration for the given
    % run/trial and builds the default report (summary slide, compare
    % slide, archived figures, requested display plots). Each step is also
    % independently callable afterward via the public methods below.
    %
    % Options (varargin, via CheckOption):
    %   'pptfn'         (default '')              passed to Cppt() when opening a new PowerPoint
    %   'Sppt'          (default 'new')           'new' | [] | existing Cppt object
    %   'run_bn'        (default 'falco_testbed_run<runnum>')
    %   'listS'         (default [])              reuse an already-loaded CRunData-subclass array
    %   'RunDataClass'  (default 'CfalcoRunData') class constructed per-iteration, via feval()
    %   'MaxIterGuess'  (default 1000)            upper probe bound for auto-discovery
    %   'excludeItnum'  (default [])              iterations excluded from itnum_best selection
    %
    % Plus any display-plot specs (bare flag names or {method, opts...}
    % cells) and their associated clim/xlim/ylim option keywords -- see
    % AddDisplayPlots.
    %
    % Example:
    %   R = CEFCReport(203, 128);
    %   R.AddDisplayPlots(110:116, 'DisplayAllInt', 'DisplayDEfields', 'DisplayCEfields');
    %   R.AddCompareSlide(124, 128);
    %   R.Save();
    %
    % See also: GenerateEFCReport_falco (backward-compatible wrapper), CRunData, CfalcoRunData, Cppt

    properties
        runnum
        TrialNum
        run_bn
        report_pn

        S               % array of CRunData-subclass objects (e.g. CfalcoRunData), all loaded iterations
        Sppt            % [] or a Cppt instance (normalized once, never a raw string)
        itnum_best      % best iteration found by most recent AddSummarySlide() call
        ExcludeItnum    % iteration numbers excluded from itnum_best selection (still shown in trend plot)
        RunDataClass    % name of class constructed per-iteration via feval() (default 'CfalcoRunData')

        probeh          % mean probe intensity per mode, from most recent AddSummarySlide() call
        rmsdDMv         % rms delta DM volts per DM, from most recent AddSummarySlide() call
        listHfig        % figure handles created by most recent AddDisplayPlots() call
    end % properties

    methods

        function R = CEFCReport(runnum, TrialNum, varargin)

            % options
            ppt_fn = CheckOption('pptfn', '', varargin{:});
            Sppt = CheckOption('Sppt', 'new', varargin{:});
            bSpptIsNew = ischar(Sppt) && strcmpi(Sppt, 'new'); % Sppt gets resolved below, save this first
            listSin = CheckOption('listS', [], varargin{:}); % if listS of CRunData-subclass for iterations already exists

            R.runnum = runnum;
            R.TrialNum = TrialNum;
            R.run_bn = CheckOption('run_bn', ['falco_testbed_run' num2str(runnum)], varargin{:});
            R.RunDataClass = CheckOption('RunDataClass', 'CfalcoRunData', varargin{:});
            R.ExcludeItnum = CheckOption('excludeItnum', [], varargin{:});

            % resolve Sppt once, correctly (bug fix: on non-PC, 'new' must
            % normalize to [], not stay a literal string -- otherwise
            % ~isempty(Sppt) checks downstream wrongly take the "add to
            % PowerPoint" branch and crash instead of falling back to
            % per-plot-type-folder image output)
            if bSpptIsNew && ispc
                R.Sppt = Cppt(ppt_fn);
            elseif bSpptIsNew
                R.Sppt = []; % non-pc: folder-per-plot-type fallback
            else
                R.Sppt = Sppt; % already [] or a Cppt object
            end

            more off

            % load all the iteration data, if not input
            % if an existing listS is passed in, this is a "reuse" call:
            % the data is already loaded and (by convention) already has a
            % summary/compare slide in the report -- just wire up R.S /
            % R.Sppt / R.report_pn and let the caller explicitly invoke
            % whichever methods it wants next (e.g. AddDisplayPlots).
            % Only a fresh auto-discovery runs the default full-report
            % sequence, so that a "later call" reusing listS/Sppt never
            % re-inserts a duplicate summary+compare slide pair.
            bReuse = ~isempty(listSin);

            if bReuse
                R.S = listSin;
            else
                R.AutoLoadIterations(varargin{:});
            end

            % path definitions
            R.report_pn = PathTranslator(fullfile(getenv("DATA_ROOT"), R.run_bn, 'reports', R.S(1).runLabel));

            if ~bReuse
                % default full-report sequence, built from public methods
                R.AddSummarySlide();
                R.AddCompareSlide();

                if bSpptIsNew
                    R.ImportArchivedFigures();
                end

                R.AddDisplayPlots(varargin{:});

                R.Save();

            end
            
            more on

        end % constructor

        function AutoLoadIterations(R, varargin)
            % AutoLoadIterations(R, varargin)
            %
            % probe itnum = 1:MaxIterGuess, constructing R.RunDataClass
            % objects via feval(), guarded by try/catch (a nonexistent
            % iteration may throw outright rather than just returning an
            % empty ImCube). Stops after max_empties successive
            % failed/empty iterations, or once the constructed object
            % reports the run's last iteration (falcoData.Itr).

            max_empties = CheckOption('max_empties', 3, varargin{:});
            MaxIterGuess = CheckOption('MaxIterGuess', 1000, varargin{:});

            mp = [];
            cnt_empty = 0;
            n_trailing_appended_empty = 0;
            Stmp = [];
            for itnum = 1:MaxIterGuess
                fprintf('reading itnum %d\n', itnum);
                try
                    Sii = feval(R.RunDataClass, R.runnum, R.TrialNum, itnum, 'mp', mp, varargin{:});
                catch ME
                    disp(ME.message);
                    cnt_empty = cnt_empty + 1;
                    if cnt_empty >= max_empties
                        break
                    end
                    continue
                end

                if isempty(Stmp)
                    Stmp = Sii;
                else
                    Stmp(end+1) = Sii; %#ok<AGROW>
                end

                if isempty(Sii.ImCube)
                    cnt_empty = cnt_empty + 1;
                    n_trailing_appended_empty = n_trailing_appended_empty + 1;
                else
                    % reset, only count successive empties
                    cnt_empty = 0;
                    n_trailing_appended_empty = 0;
                end

                if cnt_empty >= max_empties
                    if n_trailing_appended_empty > 0
                        Stmp = Stmp(1:end-n_trailing_appended_empty);
                    end
                    break
                end

                % update mp, so it is used for next iteration and avoid
                % reading config every iteration
                mp = Sii.mp;

                % out.Itr from "_snippet.mat" is last iteration
                if itnum >= Sii.falcoData.Itr
                    break
                end
            end % for itnum

            if isempty(Stmp)
                error('CEFCReport:AutoLoadIterations', ...
                    'no valid iterations found for run %d trial %d', R.runnum, R.TrialNum);
            end

            % if last iteration is empty, remove it
            if isempty(Stmp(end).ImCube)
                Stmp = Stmp(1:end-1);
            end

            R.S = Stmp;

        end % AutoLoadIterations

        function Sout = SelectIterations(R, itnumRange)
            % Sout = SelectIterations(R, itnumRange)
            %
            % filter already-loaded R.S down to itnumRange. itnumRange may
            % be numeric itnum values, or [] / omitted to select all.

            if ~exist('itnumRange', 'var') || isempty(itnumRange)
                Sout = R.S;
                return
            end
            [~, i_S] = intersect([R.S.iter], itnumRange);
            if isempty(i_S)
                error('CEFCReport:SelectIterations', 'requested itnumRange not found in loaded iterations');
            end
            Sout = R.S(i_S);

        end % SelectIterations

        function itnum_best = AddSummarySlide(R, varargin)
            % itnum_best = AddSummarySlide(R, varargin)
            %
            % 2x2 slide: NormIntensity, Beta, Probeh, RMSdDMv vs iteration.
            % Sets R.itnum_best = lowest mean normalized intensity,
            % excluding any itnum in the 'excludeItnum' option (default
            % R.ExcludeItnum) -- excluded iterations still appear in the
            % plotted trend curve but can't be selected.

            excludeItnum = CheckOption('excludeItnum', R.ExcludeItnum, varargin{:});
            S = R.S;

            hfig = figure_mxn(2,2);
            set(hfig, 'position', [322 263 1800 1000]);
            hax(1,1) = subplot(2,2,1); [~, ~, ~, itnum_best] = R.PlotNormIntensity('hfig', hfig, 'hax', hax(1,1), 'excludeItnum', excludeItnum);
            hax(1,2) = subplot(2,2,2); R.PlotBeta('hfig', hfig, 'hax', hax(1,2));
            hax(2,1) = subplot(2,2,3); probeh = R.PlotProbeh('hfig', hfig, 'hax', hax(2,1));
            hax(2,2) = subplot(2,2,4); rmsdDMv = R.PlotRMSdDMv('hfig', hfig, 'hax', hax(2,2));

            R.AddFigureToReport(hfig, 'summary', ['summary_it' num2str(S(1).iter) '_it' num2str(S(end).iter)]);

            R.itnum_best = itnum_best;
            R.probeh = probeh;
            R.rmsdDMv = rmsdDMv;

        end % AddSummarySlide

        function sImageCompareData = AddCompareSlide(R, itnum_compare1, itnum_compare2, varargin)
            % sImageCompareData = AddCompareSlide(R, itnum_compare1, itnum_compare2, varargin)
            %
            % compare full-band (mean over subbands/modes) normalized
            % intensity images -- Total (unprobed), Modulated, Unmodulated
            % -- side by side between two iterations (default: S(end).iter-8
            % and S(end).iter-4).

            S = R.S;
            if ~exist('itnum_compare1', 'var') || isempty(itnum_compare1)
                itnum_compare1 = S(end).iter - 8;
            end
            if ~exist('itnum_compare2', 'var') || isempty(itnum_compare2)
                itnum_compare2 = S(end).iter - 4;
            end

            [hfig_cmp, ~, sImageCompareData] = R.PlotCompareNormInt(itnum_compare1, itnum_compare2, varargin{:});

            % bug fix: always save an image too, not just the .mat, when
            % there's no PowerPoint to copy the figure into
            R.AddFigureToReport(hfig_cmp, 'compare_normint', ['compare_it' num2str(itnum_compare1) '_it' num2str(itnum_compare2)]);

            if ~exist(R.report_pn, 'dir'), mkdir(R.report_pn); end
            save(fullfile(R.report_pn, ['compare_normint_it' num2str(itnum_compare1) '_it' num2str(itnum_compare2) '.mat']), 'sImageCompareData');

        end % AddCompareSlide

        function ImportArchivedFigures(R)
            % ImportArchivedFigures(R)
            %
            % import .png figures already saved by falco itself into the
            % report. Only meaningful when this call created a brand-new
            % PowerPoint (guarded by the caller, e.g. the constructor).

            if isempty(R.Sppt)
                return
            end
            S = R.S;
            figures_pn = [S(1).Rundir_pn '/figures'];
            if ~exist(PathTranslator(figures_pn), 'dir')
                return
            end
            listPng = dir(fullfile(PathTranslator(figures_pn), '*.png'));
            for ii = 1:length(listPng)
                slide = R.Sppt.NewSlide(1+ii);
                fn = fullfile(listPng(ii).folder, listPng(ii).name);
                invoke(slide.Shapes, 'AddPicture', fn, true, true, 100, 100);
            end

        end % ImportArchivedFigures

        function listHfig = AddDisplayPlots(R, varargin)
            % listHfig = AddDisplayPlots(R, varargin)
            % listHfig = AddDisplayPlots(R, itnumRange, varargin)
            %
            % varargin entries can be:
            %   bare flag name (char), e.g. 'DisplayAllInt', 'DisplayCEfields', ...
            %      -> translated to {method, opt1, val1, ...} using the same
            %         clim/xlim/ylim defaults GenReport.m used to hardcode
            %   {method, opt1, val1, ...} cell, passed through as-is
            %      (today's direct GenerateEFCReport_falco calling convention)
            %
            % An optional leading numeric itnumRange restricts which
            % already-loaded R.S entries the plots run against (a pure
            % filter -- R.S itself always holds every loaded iteration).
            %
            % Options (apply only to the bare-flag-name translation):
            %   'clim_defields', []
            %   'clim_intensity', [-9 -6]
            %   'clim_probecube', 'auto'
            %   'clim_deltaDMv', []
            %   'xlim_display', -10 + 10*[-1 1]
            %   'ylim_display', 0.8 + 10*[-1 1]

            itnumRange = [];
            if ~isempty(varargin) && isnumeric(varargin{1})
                itnumRange = varargin{1};
                varargin = varargin(2:end);
            end

            Ssel = R.SelectIterations(itnumRange);

            listSpec = R.ResolveDisplaySpecs(varargin{:});

            listHfig = {};
            for iplot = 1:length(listSpec)
                listHfig{end+1} = R.CreatePlots(Ssel, listSpec{iplot}{1}, listSpec{iplot}{2:end}); %#ok<AGROW>
            end

            R.listHfig = listHfig;

        end % AddDisplayPlots

        function Save(R, varargin)
            % Save(R, varargin)
            %
            % PPTX SaveAs, with an overwrite-confirmation prompt if the
            % target file already exists. No-op if there's no open Sppt.

            if isempty(R.Sppt)
                return
            end
            fn = CheckOption('fn', PathTranslator(fullfile(getenv("DATA_ROOT"), R.run_bn, 'reports', ...
                [R.S(1).runLabel '_it' num2str(R.S(1).iter) '_' num2str(R.S(end).iter) '.pptx'])), varargin{:});

            if exist(fn, 'file')
                choice = questdlg(['File already exists: ' fn], ...
                    'File Exists', ...
                    'Replace', 'Save with new name', 'Cancel', 'Cancel');

                switch choice
                    case 'Replace'
                        delete(fn);
                        R.Sppt.Presentation.SaveAs(fn);
                        fprintf('Replaced existing file: %s\n', fn);

                    case 'Save with new name'
                        [fpath, fname, fext] = fileparts(fn);
                        datetime_str = datestr(now, 'yyyymmdd_HHMMSS');
                        fn_new = fullfile(fpath, [fname '_' datetime_str fext]);
                        R.Sppt.Presentation.SaveAs(fn_new);
                        fprintf('Saved as: %s\n', fn_new);

                    otherwise
                        fprintf('Save cancelled. File not saved.\n');
                end
            else
                R.Sppt.Presentation.SaveAs(fn);
                fprintf('Saved: %s\n', fn);
            end

        end % Save

        % Plot* helpers are public (not private) so that external callers
        % can hold a function handle to them, e.g. the backward-compat
        % GenerateEFCReport_falco.m wrapper's sOut.fPlotX fields -- MATLAB
        % checks method access at handle-creation time, so a handle to a
        % private method can't be created from outside the class.

        function [probeh, hfig, hax] = PlotProbeh(R, varargin)
            % [probeh, hfig, hax] = PlotProbeh(R, varargin)

            hax = CheckOption('hax', [], varargin{:});
            S = R.S;

            if ~isempty(hax)
                [~, ~, itnum_texp, texp] = R.PlotTexp('nodisplay', true);
                hfig = hax.Parent;
            else
                [hfig, hax, itnum_texp, texp] = R.PlotTexp();
            end

            itnum = zeros(size(S));
            probeh = cell(1, S(1).NofW);
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
            strStar = {'On-axis', 'Off-axis'};
            symbols = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', '+', 'x'};
            for imode = 1:length(probeh)
                ibnd = rem(imode-1, S(1).Nlamcorr) + 1;
                istar = floor((imode - ibnd)./S(1).Nlamcorr) + 1;
                hl(imode) = semilogy(itnum, probeh{imode}, ['-' symbols{imode}]); %#ok<NASGU>
                grid on, hold on
                legstr{imode} = ['Band ' num2str(ibnd) ', ' strStar{istar} ' Star'];
            end
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

        function [hfig, hax, itnum, texp] = PlotTexp(R, varargin)

            nodisplay = CheckOption('nodisplay', false, varargin{:}); % in case you only want the texp data
            S = R.S;

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

        function [betaused, betamin, hfig, hax] = PlotBeta(R, varargin)
            % [betaused, betamin, hfig, hax] = PlotBeta(R, varargin)
            %
            % falco does not record betamin

            hfig = CheckOption('hfig', [], varargin{:});
            hax = CheckOption('hax', [], varargin{:});
            listS = R.S;
            itnum = CheckOption('itnum', [listS.iter], varargin{:});

            % look for falcoData.ctrl.log10regHist and falcoData.ctrl.dmfacHist
            betaused = listS(end).falcoData.log10regHist(itnum);
            betamin = [];
            dmfac_used = listS(end).falcoData.ctrl.dmfacHist(itnum);

            if isempty(hfig)
                hfig = figure;
            else
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

        function [hfig, hax, han, itnum_best] = PlotNormIntensity(R, varargin)
            % [hfig, hax, han, itnum_best] = PlotNormIntensity(R, varargin)

            listS = R.S;
            hfig = CheckOption('hfig', [], varargin{:});
            hax = CheckOption('hax', [], varargin{:});
            itnum = CheckOption('itnum', [listS.iter], varargin{:});
            ylim = CheckOption('ylim', [], varargin{:});
            excludeItnum = CheckOption('excludeItnum', R.ExcludeItnum, varargin{:});

            itnum = itnum(:); % force column vector

            NInt_co = listS(1).falcoData.normIntModScore(itnum, :); % column per band*star = listS(1).NofW
            NInt_inco = listS(1).falcoData.normIntUnmodScore(itnum, :); % column per band*star = listS(1).NofW
            NInt_total = listS(1).falcoData.normIntMeasScore(itnum, :); % column per band*star
            if isfield(listS(1).mp, 'toggledMSWC') && listS(1).mp.toggledMSWC
                NInt_mean = sum(NInt_total, 2); % mean of subbands and stars. Not correct for toggled
            else
                NInt_mean = mean(NInt_total, 2); % mean of subbands and stars. Not correct for toggled
            end

            if isempty(hfig)
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

            strStar = {'On-axis', 'Off-axis'};
            for ibnd = 1:listS(1).NofW % check length of .lambda always equals Nsbd * Nstar ???
                if listS(1).NofW == listS(1).Nstar
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
            legend(legstr_total{:},legstr_unmod{:},legstr_mod{:})

            % best-iteration search: excluded iterations are NaN'd out
            % first (min() ignores NaN automatically), so they still show
            % in the trend curve above but can't be selected as best
            NInt_mean_search = NInt_mean;
            NInt_mean_search(ismember(itnum, excludeItnum)) = NaN;
            [itnum_best, NInt_total_min, NInt_inco_min, NInt_co_min, NInt_mean_min] = mindata(NInt_mean_search, itnum, mean(NInt_total, 2), mean(NInt_inco, 2), mean(NInt_co, 2), NInt_mean);
            han = FigureTitle(['Trial # ' num2str(listS(1).trialNum) '; Iter #' num2str(itnum_best) '; NI = ' num2str(NInt_mean_min,'%.1e') '; Mod = ' num2str(NInt_co_min,'%.1e ') '; Unmod = ' num2str(NInt_inco_min,'%.1e ')],'FontSize',12);

        end % PlotNormIntensity

        function [rmsdDMv, hfig, hax] = PlotRMSdDMv(R, varargin)
            % [rmsdDMv, hfig, hax] = PlotRMSdDMv(R, varargin)

            hfig = CheckOption('hfig', [], varargin{:});
            hax = CheckOption('hax', [], varargin{:});
            listS = R.S;
            itnum = CheckOption('itnum', [listS.iter], varargin{:});

            % plotting differential, start with itnum(2)
            itnum_plot = itnum(2:end);
            itnum_plot = itnum_plot(:);

            Ndm = length(listS(1).DMvCube);

            % select the unprobed DMv from each cube
            for ii = 1:length(listS)
                for idm = 1:Ndm
                    if ~isempty(listS(ii).DMvCube)
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
                for idm = 1:Ndm
                    iuse = abs(listDMv{ii+1,idm}) > 0 & abs(listDMv{ii,idm}) > 0;
                    rmsdDMv(ii,idm) = rms(listDMv{ii+1,idm}(iuse) - listDMv{ii,idm}(iuse));
                end
            end

            if isempty(hfig)
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
            for idm = 1:Ndm
                legstr{idm} = ['DM ' num2str(idm)];
            end
            legend(legstr{:}, 'Mean')
            hax = gca;

        end % PlotRMSdDMv

        function [hfig, hax, sImageData] = PlotCompareNormInt(R, itnum1, itnum2, varargin)
            % [hfig, hax, sImageData] = PlotCompareNormInt(R, itnum1, itnum2, varargin)
            %
            % compare full-band (mean over subbands/modes) normalized
            % intensity images between two iterations: 3 rows x 2 cols
            %    row 1 = Total (unprobed)
            %    row 2 = Modulated
            %    row 3 = Unmodulated
            %    col 1 = itnum1, col 2 = itnum2
            %
            % sImageData = struct with the images used in the plot, for saving to .mat
            %
            % Options:
            %     CheckOption('bLog', true, varargin{:});
            %     CheckOption('clim', [-9 -6.5], varargin{:});
            %     CheckOption('xlim', S(1).XlimDefault, varargin{:});
            %     CheckOption('ylim', S(1).YlimDefault, varargin{:});
            %     CheckOption('title_left', ['Iter #' num2str(itnum1)])
            %     CheckOption('title_right', ['Iter #' num2str(itnum2)])
            %     CheckOption('save_fn', [], varargin{:});

            S = R.S;

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

            R.AddFigureToReport(hfig, 'PlotCompareNormInt', ['Iter_' num2str(itnum1) '_iter_' num2str(itnum2)]);
        end % PlotCompareNormInt

    end % methods (public, plot helpers)

    methods (Access = private)

        function AddFigureToReport(R, hfig, category, basename)
            % AddFigureToReport(R, hfig, category, basename)
            %
            % if R.Sppt exists, append hfig as a new slide; else save
            % hfig as report_pn/category/basename.jpg and .fig.
            % Always appends (Sppt.NewSlide([]) internally, via
            % CopyFigNewSlide) -- never a hardcoded slide index.

            if isempty(hfig)
                return
            end
            if ~isempty(R.Sppt)
                R.Sppt.CopyFigNewSlide(hfig);
            else
                R.fSaveas(hfig, category, basename);
            end

        end % AddFigureToReport

        function fSaveas(R, hfig, category, basename)
            % fSaveas(R, hfig, category, basename)
            %
            % save figure to report_pn/category/basename.jpg and .fig

            fn = fullfile(R.report_pn, category, [basename '.jpg']);
            fnfig = fullfile(R.report_pn, category, [basename '.fig']);

            pn = fileparts(fn);
            if ~exist(pn, 'dir')
                mkdir(pn);
                if isunix
                    system(['chmod 775 ' pn]);
                end
            end
            saveas(hfig, fn);
            saveas(hfig, fnfig);

        end % fSaveas

        function listSpec = ResolveDisplaySpecs(R, varargin) %#ok<INUSL>
            % listSpec = ResolveDisplaySpecs(R, varargin)
            %
            % translate bare flag names to {method, opts...} cells
            % (GenReport.m's hardcoded defaults), pass {method,opts...}
            % cell-array specs through unchanged.

            clim_defields = CheckOption('clim_defields', [], varargin{:});
            clim_intensity = CheckOption('clim_intensity', [-9 -6], varargin{:});
            clim_probecube = CheckOption('clim_probecube', 'auto', varargin{:});
            xlim_display  = CheckOption('xlim_display', -10 + 10*[-1 1], varargin{:});
            ylim_display = CheckOption('ylim_display',  0.8 + 10*[-1 1], varargin{:});
            clim_deltaDMv = CheckOption('clim_deltaDMv', [], varargin{:});

            listSpec = {};
            for ii = 1:length(varargin)
                arg = varargin{ii};
                if iscell(arg)
                    listSpec{end+1} = arg; %#ok<AGROW>
                    continue
                end
                if ~ischar(arg)
                    continue % option value already consumed by CheckOption above
                end

                switch arg
                    case 'DisplayAllInt'
                        listSpec{end+1} = {'DisplayAllInt', 'clim', clim_intensity, 'DisplayRadialIntensity', false, 'xlim', xlim_display, 'ylim', ylim_display}; %#ok<AGROW>
                    case 'DisplayCEfields'
                        listSpec{end+1} = {'DisplayCEfields', 'nodisplay', false, 'PSF_thresh_nsig', 3, 'xlim', xlim_display, 'ylim', ylim_display}; %#ok<AGROW>
                    case 'DisplayDEfields'
                        listSpec{end+1} = {'DisplayDEfields', 'clim', clim_defields, 'xlim', xlim_display, 'ylim', ylim_display}; %#ok<AGROW>
                    case 'DisplayProbeCube'
                        listSpec{end+1} = {'DisplayProbeCube', 'xlim', xlim_display, 'ylim', ylim_display, 'clim', clim_probecube, 'iwv', 1}; %#ok<AGROW>
                        listSpec{end+1} = {'DisplayProbeCube', 'xlim', xlim_display, 'ylim', ylim_display, 'clim', clim_probecube, 'iwv', 2}; %#ok<AGROW>
                    case 'DisplayDMv'
                        listSpec{end+1} = {'DisplayDMv', 'climdelta', clim_deltaDMv}; %#ok<AGROW>
                    case 'DisplayDMvProbe'
                        listSpec{end+1} = {'DisplayDMvProbe'}; %#ok<AGROW>
                    case 'DisplayEfields'
                        listSpec{end+1} = {'DisplayEfields', 'clim', 1e-2*[-1 1]}; %#ok<AGROW>
                    otherwise
                        % one of the option keywords consumed above (e.g.
                        % 'clim_intensity') -- silently skip, matching the
                        % CheckOption convention used throughout this codebase
                end
            end

        end % ResolveDisplaySpecs

        function [hfig, hax, sCmetrics] = CreatePlots(R, S, sDisplayFun, varargin)
            % [hfig, hax, sCmetrics] = CreatePlots(R, S, sDisplayFun, varargin)
            %
            % create the plots for one display method across the
            % iterations in S. Some plots are differential (need a
            % reference iteration); some also produce a summary trend
            % slide/figure.

            figheight = CheckOption('figheight', 700, varargin{:}); % for ppt display
            trialname = CheckOption('trialname', '', varargin{:});

            if ~ispc && ~exist(R.report_pn, 'dir')
                mkdir(R.report_pn);
                if isunix
                    system(['chmod 775 ' R.report_pn]);
                end
            end

            % CRunData methods where the first argument is a reference iteration
            listDiff = {
                'DisplayDEfields'
                'DisplayDMv'
                };

            N = length(S);

            % if only 1 iteration, can't do displays that use differences
            if N <= 1 && any(strcmp(sDisplayFun, [listDiff(:); {'DisplayCEfields'}]))
                hfig = []; hax = []; sCmetrics = struct;
                return
            end

            hfig = [];
            switch sDisplayFun

                case 'DisplayDEfields'
                    for ii = 1:N-1
                        sMtmp = [];
                        try
                            [hfig, hax, sMtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:}, 'hfig', hfig);
                        catch ME
                            hfig = [];
                            disp(ME.message);
                            for ierr = 1:length(ME.stack)
                                disp(ME.stack(ierr));
                            end
                        end
                        if ~isempty(sMtmp)
                            sCmetrics(ii) = sMtmp;
                        end

                        if ~isempty(hfig)
                            figscale = R.CalcFigscale(hfig, figheight);
                            set(hfig, 'Position', figscale*get(hfig,'position'));

                            % if using ImageCube to toggle model/measure, make the gif
                            if ~isempty(get(hfig, 'KeyPressFcn'))
                                fungif = get(hfig, 'KeyPressFcn');

                                gif_fn = fullfile(R.report_pn, sDisplayFun, ['it_' num2str(S(ii).iter) '.gif']);
                                pn = fileparts(gif_fn);
                                if ~exist(pn, 'dir')
                                    mkdir(pn);
                                    if isunix
                                        system(['chmod 775 ' pn]);
                                    end
                                end

                                % this sends an event to hfig to create the gif
                                thisevent = struct('Modifier', 'shift', 'Key', 'g', 'gif_fn', gif_fn);
                                fungif(hfig, thisevent);

                                % check and insert gif to PowerPoint
                                if exist(gif_fn, 'file') && ~isempty(R.Sppt)
                                    R.Sppt.AddPictureNewSlide(gif_fn);
                                end

                            else
                                R.fSaveas(hfig, sDisplayFun, ['it_' num2str(S(ii).iter)]);
                            end % if ImageCube
                        end % if hfig
                    end % for ii iter

                    % plot a summary rms dE vs iteration
                    if any(strcmp({sCmetrics.type}, 'dEfields'))
                        nw = S(1).NofW; % for convenience
                        hfig_de = figure;
                        hh = semilogy([S(2:N).iter], [sCmetrics.rmsdE_t].^2, '-', ...
                            [S(2:N).iter], mean([sCmetrics.rmsdE_t].^2, 1), '-', ...
                            [S(2:N).iter], [sCmetrics.rmsdE_m].^2, '--', ...
                            [S(2:N).iter], mean([sCmetrics.rmsdE_m].^2, 1), '--');

                        grid on
                        set(hh(nw+1),'linewidth',2)
                        set(hh(end), 'linewidth',2)
                        wvstrtmp = num2str(S(1).NKTcenter(1:nw)'/S(1).NM);
                        tbwvstrtmp = char(strcat('TB', {' '}, wvstrtmp, 'nm'));
                        mowvstrtmp = char(strcat('Model', {' '}, wvstrtmp, 'nm'));
                        legend(char(tbwvstrtmp, 'TB-mean', mowvstrtmp, 'Model-Mean'))
                        xlabel('Iteration #')
                        ylabel('mean |\DeltaE|^2')
                        title(trialname, 'fontsize', 14)

                        R.AddFigureToReport(hfig_de, 'summary', ['magdE_it' num2str(S(1).iter) '_it' num2str(S(end).iter)]);
                    end

                case 'DisplayDMv'
                    for ii = 1:N-1
                        sMtmp = [];
                        try
                            [hfig, hax, sMtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:}, 'hfig', hfig);
                        catch ME
                            hfig = [];
                            disp(ME.message);
                            for ierr = 1:length(ME.stack)
                                disp(ME.stack(ierr));
                            end
                        end
                        if ~isempty(sMtmp)
                            sCmetrics(ii) = sMtmp;
                        end

                        if ~isempty(hfig)
                            figscale = R.CalcFigscale(hfig, figheight);
                            set(hfig, 'Position', figscale*get(hfig,'position'));
                            R.AddFigureToReport(hfig, sDisplayFun, ['it_' num2str(S(ii).iter)]);
                        end
                    end % for ii iter

                case 'DisplayCEfields'
                    for ii = 1:N-1
                        sCtmp = [];
                        try
                            [hfig, hax, sCtmp] = S(ii+1).(sDisplayFun)(S(ii), varargin{:}, 'hfig', hfig);
                        catch ME
                            hfig = [];
                            disp(ME.message);
                            for ierr = 1:length(ME.stack)
                                disp(ME.stack(ierr));
                            end
                        end

                        if ~isempty(sCtmp)
                            sCmetrics(ii) = sCtmp;
                        end

                        if isempty(hfig)
                            continue
                        end

                        figscale = R.CalcFigscale(hfig, figheight);
                        set(hfig, 'Position', figscale*get(hfig,'position'));
                        R.AddFigureToReport(hfig, sDisplayFun, ['it_' num2str(S(ii).iter)]);
                    end % for ii iter

                    hfig_ce = figure;
                    plotampphase([S(2:N).iter], [sCmetrics.CC], ...
                        'xlabel','Iteration #','title',[trialname ', \DeltaE Testbed Model Correlation (CC)']);

                    R.AddFigureToReport(hfig_ce, 'summary', ['CE_it' num2str(S(1).iter) '_it' num2str(S(end).iter)]);

                otherwise % one call per iteration
                    sCmetrics = struct;
                    for ii = 1:N
                        try
                            hfig = S(ii).(sDisplayFun)(varargin{:}, 'hfig', hfig);
                        catch ME
                            hfig = [];
                            disp(ME.message);
                            for ierr = 1:length(ME.stack)
                                disp(ME.stack(ierr));
                            end
                        end

                        if isempty(hfig)
                            continue
                        end

                        figscale = R.CalcFigscale(hfig, figheight);
                        set(hfig, 'Position', figscale*get(hfig,'position'));
                        R.AddFigureToReport(hfig, sDisplayFun, ['it_' num2str(S(ii).iter)]);
                    end

            end % switch

        end % CreatePlots

        function figscale = CalcFigscale(R, hfig, figheight) %#ok<INUSL>
            pos = get(hfig,'Position');
            ysize = pos(end);
            figscale = figheight/ysize;
        end % CalcFigscale

    end % methods (private)

end % classdef
