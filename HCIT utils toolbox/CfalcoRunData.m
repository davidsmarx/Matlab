classdef CfalcoRunData < CRunData
    
    properties
        % seriesNum = CRunData.runnum
        trialNum
        % iter = CRunData.iter
        
        % falco data
        runLabel;
        
        % falco data from snippet:
        falcoData
        
    end % properties
    
    
    methods
        
        function S = CfalcoRunData(seriesNum, trialNum, iter, varargin)
            
            %
            S.runnum = seriesNum;
            S.trialNum = trialNum;
            S.iter = iter;
            
            % paths to data
            S.runLabel = ['Series',num2str(seriesNum,'%04d'),'_Trial',num2str(trialNum,'%04d')];
            S.Rundir_pn = ['/home/hcit/OMC/OMC_MSWC/falco_testbed_run' num2str(seriesNum) '/data/' S.runLabel]; % for snippet file
            S.Reduced_pn = [S.Rundir_pn '/' S.runLabel];
            
            % config and snippet, S.falcoData
            S.loadRunData;
            
            % Load the FALCO probing data
            S.ReadProbeCube;
            if isempty(S.ImCube)
                % no data, empty iteration
                return
            end

            % need:
            %mp.Fend.corr.maskBool;
            %mp.Fend.score.maskBool;
            
            % Load FITS image data IncInt, IncIntEst, CohInt, ImCubeUnProb
            S.ReadImageCube;
            
            % DM cube
            S.ReadDMvCube('dm1');
            S.ReadDMvCube('dm2');
            
            % Get fits keywords
            finfo = fitsinfo(PathTranslator([S.Reduced_pn '/normI_it' num2str(S.iter) '.fits']));
            S.ReducedKeys = finfo.PrimaryData.Keywords;
            iwv = 1;
            S.NKTlower(iwv) = FitsGetKeywordVal(S.ReducedKeys,'NKTLOWER')*S.NM;
            S.NKTupper(iwv) = FitsGetKeywordVal(S.ReducedKeys,'NKTUPPER')*S.NM;
            S.NKTcenter(iwv) = mean([S.NKTlower(iwv) S.NKTupper(iwv)]);
            S.lambda = S.NKTcenter;

        end % init
        
        function loadRunData(S)
            %[falcoData,fitsData] = loadTBdata(falcoOutDir,seriesNum,trialNum_list)
            %   Function to load and unpack FALCO data from a sequence of trials
            %
            % 2022-03-31
            %   revising to use CRunData
            
            % fullfile('Y:\ln_hcit_omc\OMC_MSWC\falco_testbed_run200\data', bn, bn, ['normI_Esens_it' num2str(S.iter) '.fits']));
            
            % Load the configuration file => mp
            %config_fn = [S.Rundir_pn '/' S.runLabel '_config.m'];
            config_fn = [S.Rundir_pn '/falco_omc_config_' S.runLabel '.m'];
            config_fn = PathTranslator(config_fn);
            if ~exist(config_fn, 'file'),
                warning(['cannot find config file, using default values, ' config_fn]);
                S.ppl0 = 5.54*520/575;
                S.Nppair = 3;
                S.NofW = 1;
                S.XYlimDefault = 22;
                S.Nlamcorr = 1;
                S.XYlimDefault = 22;
                S.RminSc = 3; % = mp.Fend.score.Rin
                S.RmaxSc = 9; % = mp.Fend.score.Rout
                S.ThminSc = []; % derive from mp.Fend.score.ang & mp.Fend.sides
                S.ThmaxSc = []; % derive from mp.Fend.score.ang & mp.Fend.sides
                S.YminSc = -inf;
                S.YmaxSc = inf;
                S.XminSc = -inf;
                S.XmaxSc = inf;

            else
                %copyfile(config_fn, './config_tmp.m');
                eval('config_tmp');
                %         falcoData.Nsbp(trialIndex) = mp.Nsbp;
                %         falcoData.fracBW(trialIndex) = mp.fracBW;
                %         falcoData.lambda0(trialIndex) = mp.lambda0;
                %         falcoData.sbp_centers = mp.sbp_centers;
                %         falcoData.Nitr(trialIndex) = mp.Nitr;
                %         falcoData.si_ref(trialIndex) = mp.si_ref;
                %         falcoData.Fend = mp.Fend;
                %         falcoData.thput_metric = mp.thput_metric;
                %         falcoData.thput_radius = mp.thput_radius;
                
                S.ppl0 = mp.Fend.res;
                S.Nppair = mp.est.probe.Npairs;
                S.NofW = mp.Nsbp;
                S.Nlamcorr = mp.Nsbp;
                S.XYlimDefault = 22;
                S.RminSc = mp.Fend.score.Rin;
                S.RmaxSc = mp.Fend.score.Rout;
                S.ThminSc = []; % derive from mp.Fend.score.ang & mp.Fend.sides
                S.ThmaxSc = []; % derive from mp.Fend.score.ang & mp.Fend.sides
                S.YminSc = -inf;
                S.YmaxSc = inf;
                S.XminSc = -inf;
                S.XmaxSc = inf;
                
                S.Ndm = length(mp.dm_ind);
                
            end % if load confi
            
            % Load the "snippet" file -- struct out
            snippet_fn = [S.Rundir_pn '/' S.runLabel,'_snippet.mat'];
            S.Rundir_fn = PathTranslator(snippet_fn);
            load(PathTranslator(snippet_fn)); % out
            S.falcoData = out;
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
                        
            % need:
            
        end % end loadTBdata
        
        function ReadProbeCube(S, varargin)
            % loadProbeData(S)
            RminSc = CheckOption('RminSc', S.RminSc, varargin{:});
            RmaxSc = CheckOption('RmaxSc', S.RmaxSc, varargin{:});
            TminSc = CheckOption('TminSc', S.ThminSc, varargin{:});    
            TmaxSc = CheckOption('TmaxSc', S.ThmaxSc, varargin{:});                
            YminSc = CheckOption('YminSc', S.YminSc, varargin{:});
            YmaxSc = CheckOption('YmaxSc', S.YmaxSc, varargin{:});
            XminSc = CheckOption('XminSc', S.XminSc, varargin{:});
            XmaxSc = CheckOption('XmaxSc', S.XmaxSc, varargin{:});

            % ev
            fn = [S.Reduced_pn '/probing_data_' num2str(S.iter) '.mat'];
            if exist(PathTranslator(fn),'file')
                load(PathTranslator(fn), 'ev');
                %       imageArray: [134×134×7 double]
                %               dm1: [1×1 struct]
                %              Eest: [9428×1 double]
                %          IincoEst: [9428×1 double]
                %       IprobedMean: 4.6198e-05
                %                Im: [134×134 double]
                %             score: [1×1 struct]
                %              corr: [1×1 struct]
                %     InormProbeMax: 1.0000e-04
                %             iStar: 1
                %         ampSqMean: 4.2476e-06
                %           ampNorm: [9428×3 double]
                %        InormProbe: 5.0000e-05
                %          maskBool: [120×120 logical]
                %         amp_model: [7240×3 double]

            else
                warning('cannot find probe data');
                return
            end
            % put into CRunData object
            
            % control region and score region masks
            S.bMask = ev.maskBool; % =?= mp.Fend.corr.maskBool;

            [x, y, X, Y, R, T] = CreateGrid(S.bMask, 1./S.ppl0);
            S.bMaskSc = S.bMask & (R >= RminSc & R <= RmaxSc & Y >= YminSc & Y <= YmaxSc & X >= XminSc & X <= XmaxSc);

            % S.ProbeModel{iwl, ip} = squeeze(...
            %                         ProbeData(:,:,0*S.Nlamcorr*S.Nppair+(ip-1)*S.Nlamcorr+iwl) ...
            %                         .* exp(1i*...
            %                         ProbeData(:,:,1*S.Nlamcorr*S.Nppair+(ip-1)*S.Nlamcorr+iwl)) ...
            %                         );
            %
            %                     S.ProbeMeasAmp{iwl, ip} = real(sqrt(squeeze( ProbeData(:,:,2*S.Nlamcorr*S.Nppair + (ip-1)*S.Nlamcorr+ iwl) )));
            
            for ip = 1:S.Nppair,
                S.ProbeModel{1,ip} = zeros(size(ev.maskBool));
                S.ProbeModel{1,ip}(ev.maskBool) = ev.amp_model(:,ip); % amplitude
                S.ProbeMeasAmp{1,ip} = zeros(size(ev.maskBool));
                S.ProbeMeasAmp{1,ip}(ev.maskBool) = ev.ampNorm(:,ip)*sqrt(ev.InormProbe);
            end
            
            % probe data ev includes normalized image data
            % image cube (:,:,2*Nppair+1)
            S.ImCube = ev.imageArray;
            S.imgindex = 1; % slice of ImCube of unprobed image for each subband
            
            % E-fields, S.E_t(iwl,:,:)
            S.E_t(1,:,:) = zeros(size(S.bMask));
            S.E_t(1,S.bMask) = ev.Eest;
            
            S.E_m(1,:,:) = zeros(size(S.bMask));
            S.E_m(1,S.bMask) = ev.Esim;
            
        end % loadProbeData
        
        function ReadImageCube(S)
            % S = ReadImageCube(S)
            
            % Load FITS image data
            if exist(PathTranslator([S.Reduced_pn '/normI_Esens_it' num2str(S.iter) '.fits']),'file')
                Icoh = fitsread(PathTranslator([S.Reduced_pn '/normI_Esens_it' num2str(S.iter) '.fits']));
                Iinc = fitsread(PathTranslator([S.Reduced_pn '/normI_inco_it' num2str(S.iter) '.fits']));
                Iunpr = fitsread(PathTranslator([S.Reduced_pn '/normI_it' num2str(S.iter) '.fits']));
            else
                warning('cannot find normI image cube');
                return
            end
            
            % for each wavelength subband
            iwl = 1;
            S.IncInt{iwl} = Iinc;
            S.IncIntEst{iwl} = Iinc; S.IncIntEst{iwl}(Iinc < 0) = eps; % IncInt(IncInt < 0) = eps
            S.CohInt{iwl} = Icoh;
            S.ImCubeUnProb{iwl} = Iunpr;

            % where inc int < 0, make coh int = un probed (i.e. the whole thing            
            S.CohIntEst{iwl} = S.CohInt{iwl};
            S.CohIntEst{iwl}(S.IncInt{iwl} <= 0) = S.ImCubeUnProb{iwl}(S.IncInt{iwl} <= 0);
            
            % full band = average of subbands
            S.ImCubeUnProbFullBand = mean(cat(3, S.ImCubeUnProb{:}), 3);
            S.CohIntFullBand = mean(cat(3, S.CohIntEst{:}), 3);
            S.IncIntFullBand = mean(cat(3, S.IncInt{:}), 3);
            
        end % ReadImageCube
        
        function ReadDMvCube(S, whichdm)
            % ReadDMvCube(S, whichdm)
            % whichdm = 'dm1', or 'dm2'
                        
            %             %
            %             flnm = PathTranslator([S.Reduced_pn '/' whichdm '_model_it',num2str(S.iter),'.fits']);
            %             if exist(flnm, 'file')
            %                 dmSurf = fitsread(flnm);
            %             else
            %                 S.DMvCube{end+1} = zeros(48);
            %             end
            
            flnm = PathTranslator([S.Reduced_pn '/' whichdm '_Vbias.fits']);
            if exist(flnm,'file')
                DMdata.dmVbias = fitsread(flnm);
                flnm = PathTranslator([S.Reduced_pn '/' whichdm '_V_it',num2str(S.iter),'.fits']);
                dmV = fitsread(flnm);
                S.DMvCube{end+1} = DMdata.dmVbias + dmV;
            else
                S.DMvCube{end+1} = zeros(48);
            end
            
        end % ReadDMvCube
        
    end % methods
    
end % classdef