classdef CfalcoRunData < CRunData
    
    properties
        % seriesNum = CRunData.runnum
        trialNum
        % iter = CRunData.iter
        listValidIter % list of valid iterations from this trialNum
        
        % falco data
        runLabel;
        
        % falco data from snippet:
        falcoData
        
        % mp
        mp

        % multi star parameters
        Nstar
        Nmodes % = Nstar * Nsbnd

    end % properties
    
    % list of things to fix after next trial with changes in falco:
    % loadRunData : get mp from .mat file
    % line 133 : get Nstar from mp
    % line 321 : Esim for each iMode
    
    methods
        
        function S = CfalcoRunData(seriesNum, trialNum, iter, varargin)
            %
            
            % options
            mp = CheckOption('mp', [], varargin{:});
            run_pn = CheckOption('run_pn', 'falco_testbed_run', varargin{:});
            
            %
            S.runnum = seriesNum;
            S.trialNum = trialNum;
            S.iter = iter;
            
            % paths to data
            switch S.runnum
                case {200, 201, 202}
                    S.runLabel = ['Series',num2str(seriesNum,'%04d'),'_Trial',num2str(trialNum,'%04d')];
                    S.Rundir_pn = ['/home/hcit/OMC/OMC_MSWC/' run_pn num2str(seriesNum) '/data/' S.runLabel]; % for snippet file
                    S.Reduced_pn = [S.Rundir_pn '/' S.runLabel];
                
                case 250
                    S.runLabel = ['Series',num2str(seriesNum,'%04d'),'_Trial',num2str(trialNum,'%04d')];
                    S.Rundir_pn = PathTranslator(['/proj/mcb/data/EPIC/' run_pn num2str(seriesNum) '/data/' S.runLabel]); % for snippet file
                    S.Reduced_pn = [S.Rundir_pn '/' S.runLabel];

                case 251 % OMC_EPIC model
                    S.runLabel = ['Series',num2str(seriesNum,'%04d'),'_Trial',num2str(trialNum,'%04d')];
                    S.Rundir_pn = PathTranslator(['/home/dmarx/links/HCIT/OMC_EPIC_MODEL_DATA/' run_pn num2str(seriesNum) '/data/' S.runLabel]); % for snippet file
                    S.Reduced_pn = [S.Rundir_pn '/' S.runLabel];
                    
                otherwise
                    error(['unknown seriesNum ' num2str(S.runnum)]);
            end

            % if given mp, perhaps read in another iteration of this trial
            if ~isempty(mp)
                S.mp = mp;
            end
            
            % config and snippet, S.falcoData
            S.loadRunData;
            
            % Load the FALCO probing data
            S.ReadProbeCube;
            if isempty(S.ImCube)
                % no data, empty iteration
                return
            end

            
            % Load FITS image data IncInt, IncIntEst, CohInt, ImCubeUnProb
            S.ReadImageCube;
            
            % DM cube
            S.ReadDMvCube('dm1');
            S.ReadDMvCube('dm2');
            
            % Get fits keywords
            finfo = fitsinfo(PathTranslator([S.Reduced_pn '/normI_it' num2str(S.iter) '.fits']));
            S.ReducedKeys = finfo.PrimaryData.Keywords;

            % fits headers only contain NKT values at the time the fits is
            % created, get subband centers from mp
            for iwv = 1:S.NofW                
                S.NKTcenter(iwv) = S.mp.sbp_centers(iwv);
            end
            S.lambda = S.NKTcenter;

        end % init
        
        function loadRunData(S, varargin)
            %[falcoData,fitsData] = loadTBdata(falcoOutDir,seriesNum,trialNum_list)
            %   Function to load and unpack FALCO data from a sequence of trials
            %
            % 2022-03-31
            %   revising to use CRunData
            
            mp = CheckOption('mp', S.mp, varargin{:});
            
            % fullfile('Y:\ln_hcit_omc\OMC_MSWC\falco_testbed_run200\data', bn, bn, ['normI_Esens_it' num2str(S.iter) '.fits']));
            
            if isempty(mp),
                % Load the configuration file => mp
                %config_fn = [S.Rundir_pn '/' S.runLabel '_config.m'];
                %%%% change this to read from .mat
                config_fn = [S.Rundir_pn '/falco_omc_config_' S.runLabel '.mat'];
                config_fn = PathTranslator(config_fn);
                if ~exist(config_fn, 'file'),
                    error(['cannot find config file, using default values, ' config_fn]);
                end

                % % create mp
                % copyfile(config_fn, './config_tmp.m');
                % eval('config_tmp');
                % mp = falco_flesh_out_workspace(mp);

                % load mp, copy local then load is many times faster than load from s383 server
                copyfile(config_fn, './config_tmp.mat');
                mp = load('./config_tmp.mat');

            end % if isempty(mp)
            
            %S.Nppair = mp.est.probe.Npairs;
            S.Nstar = 1; %2; % hard code for now
            S.Nlamcorr = mp.Nsbp; % 
            S.NofW = S.Nstar * S.Nlamcorr; % use NofW as all to fool CRunData
            
            %                 S.RminSc = Fend.score.Rin;
            %                 S.RmaxSc = Fend.score.Rout;
            %                 S.ThminSc = []; % derive from Fend.score.ang & Fend.sides
            %                 S.ThmaxSc = []; % derive from Fend.score.ang & Fend.sides
            %                 S.YminSc = -inf;
            %                 S.YmaxSc = inf;
            %                 S.XminSc = -inf;
            %                 S.XmaxSc = inf;

            %             % use falco to generate ctrl and score region masks
            %             score.pixresFP = Fend.res;
            %             score.rhoInner = Fend.score.Rin;
            %             score.rhoOuter = Fend.score.Rout;
            %             score.angDeg = Fend.score.ang;
            %             score.whichSide = Fend.sides;
            %             score.shape = Fend.shape;
            %             score.xiOffset = Fend.xiOffset;
            %             score.etaOffset = Fend.etaOffset;
            %             score.FOV = Fend.FOV;
            %             [S.bMaskSc, xis, etas] = falco_gen_SW_mask(score);
            %             %figure, imageschcit(bMaskSc)
            
            %             corr.pixresFP = Fend.res;
            %             corr.rhoInner = d.Fend.corr.Rin;
            %             corr.rhoOuter = Fend.corr.Rout;
            %             corr.angDeg = Fend.corr.ang;
            %             corr.whichSide = Fend.sides;
            %             corr.shape = Fend.shape;
            %             corr.xiOffset = Fend.xiOffset;
            %             corr.etaOffset = Fend.etaOffset;
            %             corr.FOV = Fend.FOV;
            %             [S.bMask] = falco_gen_SW_mask(corr);
            %figure, imageschcit(bMaskSc)
            
            S.Ndm = length(mp.dm_ind);

            % store mp in the class instance
            S.mp = mp;
                
            % Load the "snippet" file -- struct out
            snippet_fn = [S.Rundir_pn '/' S.runLabel,'_snippet.mat'];
            S.Rundir_fn = PathTranslator(snippet_fn);
            load(PathTranslator(snippet_fn), 'out'); % out
            S.falcoData = out;
            S.listValidIter = find(out.InormHist > 0);
            % struct:
            %                  Nitr: 150
            %           log10regHist: [150×1 double]
            %                   ctrl: [1×1 struct]
            %                    dm1: [1×1 struct]
            %                    dm2: [1×1 struct]
            %                    dm8: [1×1 struct]
            %                    dm9: [1×1 struct]
            %                  Zsens: [0×0×150 double]
            %      complexProjection: [149×1 double]
            %     complexCorrelation: [149×1 double]
            %              InormHist: [151×1 double]
            %           IrawCorrHist: [151×1 double]
            %          IrawScoreHist: [151×1 double]
            %           IestCorrHist: [150×1 double]
            %          IestScoreHist: [150×1 double]
            %          IincoCorrHist: [150×1 double]
            %         IincoScoreHist: [150×1 double]
            %        normIntMeasCorr: [150×1 double]
            %       normIntMeasScore: [150×1 double]
            %         normIntModCorr: [150×1 double]
            %        normIntModScore: [150×1 double]
            %       normIntUnmodCorr: [150×1 double]
            %      normIntUnmodScore: [150×1 double]
            %                  thput: [151×1 double]
            %                   Fend: [1×1 struct]
            %          serialDateVec: [150×1 double]
            %              sbp_width: 0
            %      tb_report_initial: [1×1 struct]
            %                    Itr: 116
            %          datetimeArray: [150×1 datetime]
            %           InormHist_tb: [1×1 struct]
            %            EforSpectra: {1×116 cell}
            %              smspectra: {1×116 cell}
            %                     sm: {1×116 cell}
            %                 alpha2: {1×116 cell}
        
            % rectangle score region for drawing
            S.ppl0 = S.falcoData.Fend.res;
            xminmax = S.mp.Fend.xiOffset + S.mp.Fend.score.Rout*[-1 1];
            yminmax = S.mp.Fend.etaOffset + S.mp.Fend.score.Rout*[-1 1];
            S.YminSc = yminmax(1);
            S.YmaxSc = yminmax(2);
            S.XminSc = xminmax(1);
            S.XmaxSc = xminmax(2);

            S.bMask = S.falcoData.Fend.corr.maskBool;
            S.bMaskSc = S.falcoData.Fend.score.maskBool;
            
            % default XYlim for display dark zone
            [~, ~, Xld, Yld] = CreateGrid(S.falcoData.Fend.corr.maskBool, 1./S.falcoData.Fend.res);
            S.XlimDefault = [min(Xld(S.falcoData.Fend.corr.maskBool)) max(Xld(S.falcoData.Fend.corr.maskBool))] + 2*[-1 1];
            S.YlimDefault = [min(Yld(S.falcoData.Fend.corr.maskBool)) max(Yld(S.falcoData.Fend.corr.maskBool))] + 2*[-1 1];
            
        end % end loadTBdata
        
        function ReadProbeCube(S, varargin)
            % loadProbeData(S)
            %             RminSc = CheckOption('RminSc', S.RminSc, varargin{:});
            %             RmaxSc = CheckOption('RmaxSc', S.RmaxSc, varargin{:});
            %             TminSc = CheckOption('TminSc', S.ThminSc, varargin{:});    
            %             TmaxSc = CheckOption('TmaxSc', S.ThmaxSc, varargin{:});                
            %             YminSc = CheckOption('YminSc', S.YminSc, varargin{:});
            %             YmaxSc = CheckOption('YmaxSc', S.YmaxSc, varargin{:});
            %             XminSc = CheckOption('XminSc', S.XminSc, varargin{:});
            %             XmaxSc = CheckOption('XmaxSc', S.XmaxSc, varargin{:});

            ev = CheckOption('ev', [], varargin{:});
            
            % ev
            if isempty(ev)
                fn = [S.Reduced_pn '/probing_data_' num2str(S.iter) '.mat'];
                if exist(PathTranslator(fn),'file')
                    load(PathTranslator(fn), 'ev');
                    %       imageArray: [134×134×7xNbnds double]
                    %               dm1: [1×1 struct]
                    %              Eest: [9428×Nsbnds double] % columns = subbands and stars
                    %          IincoEst: [9428×Nsbnds double]
                    %       IprobedMean: 4.6198e-05
                    %                Im: [134×134 double]
                    %             score: [1×1 struct]
                    %              corr: [1×1 struct]
                    %     InormProbeMax: 1.0000e-04
                    %             iStar: 2
                    %         ampSqMean: [1936×2×2 double]
                    %           ampNorm: [1936×2×2 double]
                    %        InormProbe: [1936×2×2 double]
                    %         amp_model: [1936×2×2 double]
                    %          maskBool: [360×360 logical]
                    %           condnum: [1.2730 1.0788 1.0747 2.1099 1.3300 1.1550 1.0561 1.1046 1.1847 1.0527 … ] (1×1936 double)
                    %              Esim: [1936×1 double]
                    
                else
                    warning(['cannot find probe data, iter#' num2str(S.iter)]);
                    return
                end
            end
            % put into CRunData object
            
            % control region and score region masks
            % were defined in loadRunData
            
            %             [x, y, X, Y, R, T] = CreateGrid(S.bMask, 1./S.ppl0);
            %             S.bMaskSc = S.bMask & (R >= RminSc & R <= RmaxSc & Y >= YminSc & Y <= YmaxSc & X >= XminSc & X <= XmaxSc);

            % S.ProbeModel{iwl, ip} = squeeze(...
            %                         ProbeData(:,:,0*S.Nlamcorr*S.Nppair+(ip-1)*S.Nlamcorr+iwl) ...
            %                         .* exp(1i*...
            %                         ProbeData(:,:,1*S.Nlamcorr*S.Nppair+(ip-1)*S.Nlamcorr+iwl)) ...
            %                         );
            %
            %                     S.ProbeMeasAmp{iwl, ip} = real(sqrt(squeeze( ProbeData(:,:,2*S.Nlamcorr*S.Nppair + (ip-1)*S.Nlamcorr+ iwl) )));

            % check that ev mask is same as corr mask
            if ~isequal(ev.maskBool, S.bMask),
                warning('ev.maskBool is not same as corr maskBool');
            end
            
            [Npx_tmp, S.Nppair, S.Nmodes] = size(ev.amp_model);
            % check: Nmodes = S.Nstar * S.Nlamcorr = S.NofW
            for iMode = 1:S.NofW % really Nmodes
            for ip = 1:S.Nppair,
                S.ProbeModel{iMode,ip} = zeros(size(ev.maskBool));
                S.ProbeModel{iMode,ip}(ev.maskBool) = ev.amp_model(:,ip, iMode); % amplitude
                S.ProbeAmp{iMode,ip} = zeros(size(ev.maskBool));
                S.ProbeAmp{iMode,ip}(ev.maskBool) = ev.ampNorm(:,ip, iMode) .* sqrt(ev.InormProbe(:, ip, iMode));
                % until we get ProbeMeasAmp from images:
                S.ProbeMeasAmp{iMode, ip} = S.ProbeAmp{iMode, ip};
            end
            end
            
            % probe data ev includes normalized image data
            % image cube (:,:, Nlamcorr * (2*Nppair+1)) % imageArray: [360×360×5×2 double]
            [ny, nx, nPrtmp, nModetmp] = size(ev.imageArray);
            if ~isequal(nPrtmp, 2*S.Nppair+1), error('ev.imageArray size inconsistent'); end
            if ~isequal(nModetmp, S.Nmodes), error('ev.imageArray size inconsistent'); end
            S.ImCube = reshape(ev.imageArray, ny, nx, []);
            %S.imgindex = 1; % slice of ImCube of unprobed image for each subband
            
            % replaces bMask defined in loadRunData, need to reconcile
            S.bMask = ev.maskBool;
            
            % E-fields, S.E_t(iwl,:,:)
            for iMode = 1:S.Nmodes
                S.E_t(iMode,:,:) = zeros(size(S.bMask));
                S.E_t(iMode,S.bMask) = ev.Eest(:,iMode);
            
                if isfield(ev, 'Esim')
                    S.E_m(iMode,:,:) = zeros(size(S.bMask));
                    S.E_m(iMode, S.bMask) = ev.Esim(:, iMode); % iMode);
                end
        
            end % for each Mode (subband, star)


        end % loadProbeData
        
        function ReadImageCube(S)
            % S = ReadImageCube(S)
            
            % Load FITS image data
            if exist(PathTranslator([S.Reduced_pn '/normI_Esens_it' num2str(S.iter) '.fits']),'file')
                Icoh = fitsread(PathTranslator([S.Reduced_pn '/normI_Esens_it' num2str(S.iter) '.fits'])); % image cube slices = subbands
                Iinc = fitsread(PathTranslator([S.Reduced_pn '/normI_inco_it' num2str(S.iter) '.fits'])); % image cube slices = subbands
                Iunpr = fitsread(PathTranslator([S.Reduced_pn '/normI_it' num2str(S.iter) '.fits'])); % total summed intensity across subbands
            else
                warning('cannot find normI image cube');
                return
            end
            
            % for each wavelength subband
            for iMode = 1:S.Nmodes
                S.IncInt{iMode} = Iinc(:,:,iMode);
                S.CohInt{iMode} = Icoh(:,:,iMode);
                S.IncIntEst{iMode} = Iinc(:,:,iMode); S.IncIntEst{iMode}(Iinc(:,:,iMode) < 0) = eps; % IncInt(IncInt < 0) = eps
                S.CohIntEst{iMode} = Icoh(:,:,iMode);
                S.ImCubeUnProb{iMode} = S.ImCube(:, :, (iMode - 1)*(2*S.Nppair + 1) + 1);
                
                %%%%% revisit:
                % % where inc int < 0, make coh int = un probed (i.e. the whole thing
                % S.CohIntEst{iMode} = S.CohInt{iMode};
                % S.CohIntEst{iMode}(S.IncInt{iMode} <= 0) = S.ImCubeUnProb{iMode}(S.IncInt{iMode} <= 0);
            end

            % full band = average of subbands
            % = mean(cat(3, S.ImCubeUnProb{:}), 3);
            S.CohIntFullBand = mean(cat(3, S.CohIntEst{:}), 3);
            S.IncIntFullBand = mean(cat(3, S.IncInt{:}), 3);
            S.ImCubeUnProbFullBand = Iunpr;

        end % ReadImageCube
        
        function ReadDMvCube(S, whichdm)
            % ReadDMvCube(S, whichdm)
            % whichdm = 'dm1', or 'dm2'
                        
            % 
            flnm = PathTranslator([S.Reduced_pn '/' whichdm '_Vbias.fits']);
            if exist(flnm,'file')
                dmVbias = fitsread(flnm);
            else
                dmVbias = zeros(48);
            end
            
            % if this dm was used for probes, get dmv cube from probe
            if num2str(S.mp.est.probe.whichDM) == whichdm(end)
                
                % ev
                fn = [S.Reduced_pn '/probing_data_' num2str(S.iter) '.mat'];
                load(PathTranslator(fn), 'ev');
                dmV_total = dmVbias + ev.(whichdm).Vall;
                
            else
                
                flnm = PathTranslator([S.Reduced_pn '/' whichdm '_V_it',num2str(S.iter),'.fits']);
                if exist(flnm, 'file')
                    dmV = fitsread(flnm);
                else
                    dmV = zeros(48);
                end
                dmV_total = dmVbias + dmV;

            end

            % rotating, etc. should be only for display
            %             % rotate, flip according to dm registration
            %             % not sure of all the possible values of orientation
            %             switch S.mp.(whichdm).orientation
            %                 case 'flipxrot180'
            %                     dmV_total = flipud(dmV_total);
            %
            %                 otherwise
            %                     error(['unknown DM orientation: ' S.mp.(whichdm).orientation]);
            %             end
                
            % add to DMvCube
            S.DMvCube{end+1} = dmV_total;
            
        end % ReadDMvCube
        
        function [hfig, hax, sMetrics] = DisplayDMv(S, dmvref, varargin)
            cOrientation = cell([1 S.Ndm]);
            for idm = 1:S.Ndm
                whichdm = ['dm' num2str(S.mp.dm_ind(idm))];
                cOrientation{idm} = S.mp.(whichdm).orientation;
            end
            [hfig, hax, sMetrics] = DisplayDMv@CRunData(S, dmvref, varargin{:}, ...
                'applyorientation', cOrientation);
        end
        
    end % methods
    
end % classdef