classdef CRunData < handle & CConstants
% class CRunData
%
% S = CRunData(runnum, iter, sOptions)
%     default  sOptions = struct(...
%                 'Results_pn', ['../results/run' num2str(runnum, '%03d') '/'] ...
%                 ,'Rundir_pn', [sOptions.results_pn 'rundir/'] ...
%                 ,'Reduced_pn', [sOptions.results_pn 'rundir/reduced/'] ...
%                 ,'PCtemp_pn', 'C:\Users\dmarx\HCITdatatemp\' ...
%                 ,'debug',false ...
%                 );
% 
% sContrast = GetContrast(S, 'RminSc',6.3,'RmaxSc',7.3)
%    sContrast = 
%     total_lam: [5x1 double]  (Norm Int over whole control region)
%     score_lam: [5x1 double]  (Norm Int over requested RminSc to RmaxSc)
%     contr_lam: [5x1 double]  (Contrast over requested RminSc to RmaxSc)
%          mean: 1.5877e-08    (Mean Contrast RminSc to RmaxSc)
%      inco_lam: [5x1 double]  (Incoherent Contrast RminSc to RmaxSc)
%        co_lam: [5x1 double]  (Coherent Contrast RminSc to RmaxSc)
%     inco_mean: 9.0238e-09    (Mean Incoherent Contrast RminSc to RmaxSc)
%       co_mean: 1.4262e-08    (Mean Coherent Contrast RminSc to RmaxSc)
%
%
% [hfig, hax] = DisplayImCubeUnProb(S, options{:})
%                 'bLog', false
%                 'drawRadii', []
%                 'clim', []
%
% [hfig, him, Im] = DisplayImCubeImage(S, imnum)
% imnum is: S.ImCube(:,:,imnum)
%
% [hfig, hax] = DisplayEfields(S, iwvplot)
%
% [hfig, hax] = DisplayProbeCube(S, iwvplot)
%
% [hfig, hax] = DisplayIncInt(S, varargin)
% options:
%    'IncIntType', default = 'normal', 'est', 'mix'
%    'drawradii', default = []
%    'bLog', default = false
%    'clim'
%    'hax'
%
% [hfig, hax] = DisplayCohInt(S, varargin)
% options:
%    'drawradii'
%    'bLog' (default = true)
%    'clim'
%    'hax'
%
% [hfig, haxlist] = DisplayAllInt(S, varargin)
%    display large table of unprobed, coh int, inc int images
%
% [hfig, haxlist] = DisplayIncCohInt(S, varargin)
%    create large display of:
%      CohInt
%      IncInt Estimated
%      IncInt Mix (not estimated)
%
% [hfig, hax] = DisplayDMv(S)
% [hfig, hax] = DisplayDMv(S, Sref)
% [hfig, hax] = DisplayDMv(S, DM1ref_fits_fn, DM2ref_fits_fn)
%
% [hfig, ha] = DisplayDEfields(S, Sref, hfig)
%

    properties
        
        % SPC disc:
        lambda          = [542 553 565 577 588]*CConstants.NM;
        %ilamcorr        = [0, 1, 2, 3, 4];
        ppl0            

        RminSc          = 6.3; % lam/D
        RmaxSc          = 19.5;
        ThminSc         = [];
        ThmaxSc         = [];
        Nbscan          = 6;

        Results_pn = '';
        Rundir_pn  = 'rundir/';         % always relative to Results_pn
        Reduced_pn = 'rundir/reduced/'; % always relative to Results_pn
        PCtemp_pn  = 'C:\Users\dmarx\HCITdatatemp\';
        Rundir_fn  = ''; % filenames include full path assembled in the init 
        Reduced_fn = '';
        
        runnum
        iter
        
        % results:
        IncInt
        IncIntEst   % part of inc int where all probeamp > 0
        IncIntMix   % part of UnProbed Image where any probeamp <= 0

        CohInt
        E_t
        E_m
        dE_ro
        dE_r1
        bPampzero
        
        ProbeAmp % measured
        
        mdMask
        bMask
        
        ImCube
        ImCubeUnProb    % {iwv}
        ImCubeDelProb   % {iwv,ipr}
        ImCubeSigProb   % {iwv,ipr} = sigtbw  in tbif.task.probes
        ImKeys
        ReducedKeys
        NofW
        NKTupper
        NKTlower
        NKTcenter
        Nppair
        NumImProbe % # of images per wavelength if probed
        Nlamcorr
        imgindex
        Contrast   % mean contrast each wavelength
        Contrast_total
        Contrast_inco
        Contrast_co
        Sthpt = []; % struct with fields Sthpt.fovx(:), Sthpt.fovy(:), Sthpt.thpt(:)
        
        Ndm    % # of DM's, see ReadDMvCube
        DMvCube
        ProbeModel
        ProbeMeasAmp
        ProbeMeasCross
        ProbeRes
        ProbeCon
        
        rplot
        IntRad
        XYlimDefault = []; % lam/D for all imageschcit
        DrawradiiDefault = [];
        DrawthetaDefault = [];

        debug = false;
        
    end % properties
        
    % reduced fits:
    % Primary: reducedcube
    %   Image: maskcube
    %   Image: probecube
    
    methods
        function S = CRunData(runnum, iter, sOptin)

            S.runnum = runnum;
            S.iter = iter;

            % default paths based on runnum
            % 'rundir_pn' and 'reduced_pn' are always relative to
            % Results_pn
            switch S.runnum,
                case 603, % SPC_disc
                    S.Results_pn = '/home/dmarx/HCIT/SPC_disc/hcim_testbed_20170705/results/run603/';
                    S.XYlimDefault = 22;
                    
                    throughput_fn = '/home/dmarx/HCIT/SPC_disc/hcim_testbed_20170705/results/Throughput_20171003T122253.mat';
                    S.Sthpt = load(PathTranslator(throughput_fn));

                case 604,
                    S.Results_pn = '/home/dmarx/HCIT/SPC_disc/hcim_testbed_20170705/results/run604/';
                    S.XYlimDefault = 22;
                case 0
                    S.Results_pn = '/home/dmarx/ln_dst_data/EFC/HLC/run000/';
                    S.XYlimDefault = 12;
                case 606,
                    S.Results_pn = '/home/dmarx/ln_mcb_data/EFC/SPC/run606/';
                    S.ppl0 = 6.09; % MCB SPC from config_MCB_SPC_20181015.py
                    S.XYlimDefault = 12;
                    S.DrawradiiDefault = [2.6 9.0];
                    S.DrawthetaDefault = 65*[-0.5 0.5]*CConstants.P;

                    S.RminSc    = 2.6; % lam/D
                    S.RmaxSc    = 9.0;
                    S.ThminSc   = (90 - 32.5)*CConstants.P;
                    S.ThmaxSc   = (90 + 32.5)*CConstants.P;

                    pntmp = '/home/dmarx/HCIT/MCB_SPC/hcim_testbed_run606/results/FOVThroughputMap/';
                    ThptCal_fn = [pntmp 'fov_20181102T110148_Mjk_Thpt.mat'];
                    % ThptCal_fn = [pntmp 'fov_20181102T110148_Flux_Thpt.mat'];
                    %S.Sthpt.fovx(:), S.Sthpt.fovy(:), S.Sthpt.thpt(:)
                    S.Sthpt = load(PathTranslator(ThptCal_fn));
                    S.Sthpt.ThptCal_fn = ThptCal_fn;

                case 607,
                    S.Results_pn = '/home/dmarx/HCIT/MCB_SPC/hcim_model2_20181021/results/run607/';
                    S.ppl0 = 6.09;
                    S.XYlimDefault = 12;
                    S.DrawradiiDefault = [2.6 9.0];
                    S.DrawthetaDefault = 65*[-0.5 0.5]*CConstants.P;
                    
                    S.RminSc    = 2.6; % lam/D
                    S.RmaxSc    = 9.0;
                    S.ThminSc   = (90 - 32.5)*CConstants.P;
                    S.ThmaxSc   = (90 + 32.5)*CConstants.P;

                case 608,
                    S.Results_pn = '/home/dmarx/ln_mcb_data/EFC/SPC/run608/';
                    S.ppl0 = 6.13; % MCB SPC from config_MCB_SPC_20181015.py
                    S.XYlimDefault = 12;
                    S.DrawradiiDefault = [2.6 9.0];
                    S.DrawthetaDefault = 65*[-0.5 0.5]*CConstants.P;

                    S.RminSc    = 2.6; % lam/D
                    S.RmaxSc    = 9.0;
                    S.ThminSc   = (90 - 32.5)*CConstants.P;
                    S.ThmaxSc   = (90 + 32.5)*CConstants.P;

                    pntmp = '/home/dmarx/HCIT/MCB_SPC/hcim_testbed_run606/results/FOVThroughputMap/';
                    ThptCal_fn = [pntmp 'fov_20181102T110148_Mjk_Thpt.mat'];
                    % ThptCal_fn = [pntmp 'fov_20181102T110148_Flux_Thpt.mat'];
                    %S.Sthpt.fovx(:), S.Sthpt.fovy(:), S.Sthpt.thpt(:)
                    S.Sthpt = load(PathTranslator(ThptCal_fn));
                    S.Sthpt.ThptCal_fn = ThptCal_fn;

                otherwise
                    error('unrecognized runnum');
            end

            % replace default values with given options
            if ~exist('sOptin','var'), sOptin = struct(); end
            ffields = fieldnames(sOptin);
            for ii = 1:length(ffields),
                S.(ffields{ii}) = sOptin.(ffields{ii});
            end
     
            % build paths and filenames for the data
            % basename, such as 'run603it00000.fits'
            s_bn = ['run' num2str(S.runnum,'%03d') 'it' num2str(iter,'%05d') '.fits'];
            S.Reduced_fn = PathTranslator([S.Results_pn S.Reduced_pn s_bn]);
            S.Rundir_fn = PathTranslator([S.Results_pn S.Rundir_pn s_bn]);

            % 'local' paths are either the same (linux cluster)
            %         or local PC drive. If PC, we will copy data
            %         locally to save unzip and read time
            if ispc,
                sReduced_local_fn = [S.PCtemp_pn S.Reduced_pn s_bn];
                sRundir_local_fn  = [S.PCtemp_pn S.Rundir_pn s_bn];
            else
                sReduced_local_fn = S.Reduced_fn;
                sRundir_local_fn  = S.Rundir_fn;
            end
            
            % unzip the fits files if necessary
            % first check if 'local' fits files exists,
            if ~exist(sReduced_local_fn, 'file'),
                sReduced_gz = [sReduced_local_fn '.gz'];
                if ~exist(sReduced_gz, 'file') && ispc,
                    copyfile([S.Reduced_fn '.gz'], sReduced_gz);
                end
                gunzip(sReduced_gz);
                % now reduced data should be unzipped, and 
                % S.Reduced_fn = the 'local' unzipped data file
            end
            S.Reduced_fn = sReduced_local_fn;

            % repeat local logic for rundir
            if ~exist(sRundir_local_fn, 'file'),
                sRundir_gz = [sRundir_local_fn '.gz'];
                if ~exist(sRundir_gz, 'file') && ispc,
                    copyfile([S.Rundir_fn '.gz'], sRundir_gz);
                end
                gunzip(sRundir_gz);
                % now rundir data should be unzipped, and 
                % S.Rundir_fn = the 'local' unzipped data file
            end
            S.Rundir_fn = sRundir_local_fn;
               
            
            finfo = fitsinfo(S.Rundir_fn);
            S.ImKeys = finfo.PrimaryData.Keywords;
            
            if ~isempty(S.Reduced_fn),
                finfo = fitsinfo(S.Reduced_fn);
                S.ReducedKeys = finfo.PrimaryData.Keywords;
            end
            
            % number of probes, and images per wavelength if probed
            % from lyotserver_IFS.py
            S.Nppair = (FitsGetKeywordVal(S.ImKeys,'NUMIM')+FitsGetKeywordVal(S.ImKeys,'PROBEF'));
            S.NumImProbe = 2*S.Nppair+1;
            % should check S.ImKeys, 'NOPAIRx' to see if wavelength x is probed

            % Norm Intensity from first, non-probed image, for each wave
            S.NofW = FitsGetKeywordVal(S.ImKeys,'NCOLOR');
            if isempty(S.NofW), S.NofW = FitsGetKeywordVal(S.ImKeys,'NCHANNEL'); end
            if isempty(S.NofW), error('keyword NCOLOR and NCHANNEL missing'); end
            S.Nlamcorr = S.NofW; % not strictly corect, should come from reduced or config length(S.ilamcorr);
            S.imgindex = 1:S.NumImProbe:S.NofW*S.NumImProbe;
            
            % get actual wavelengths from the camera fits files
            for iwv = 1:S.NofW,
               fn = FitsGetKeywordVal(S.ImKeys, ['C' num2str(iwv-1) 'P0J0']) ;
               if isempty(fn) || ~exist(PathTranslator(fn),'file'),
                   warning(['cannot open ' fn]);
                   continue
               end
               finfo = fitsinfo(PathTranslator(fn));
               if ~isempty(FitsGetKeywordVal(finfo.PrimaryData.Keywords,'NKTLOWER')),
                   S.NKTlower(iwv) = FitsGetKeywordVal(finfo.PrimaryData.Keywords,'NKTLOWER')*S.NM;
                   S.NKTupper(iwv) = FitsGetKeywordVal(finfo.PrimaryData.Keywords,'NKTUPPER')*S.NM;
                   S.NKTcenter(iwv) = mean([S.NKTlower(iwv) S.NKTupper(iwv)]);
               end
            end
            
            S.lambda = S.NKTcenter;
            
        end % function CRunData
        
        function Contrast = GetContrast(S, varargin)
            if isempty(S.bMask),
                S.ReadMaskCube;
            end

            if isempty(S.ImCube),
                S.ReadImageCube;
            end

            % parse options
            % default values from the class
            %val = CheckOption(sOpt, valDefault, varargin)
            RminSc = CheckOption('RminSc', S.RminSc, varargin{:});
            RmaxSc = CheckOption('RmaxSc', S.RmaxSc, varargin{:});
            TminSc = CheckOption('TminSc', S.ThminSc, varargin{:});    
            TmaxSc = CheckOption('TmaxSc', S.ThmaxSc, varargin{:});                
            
            [x, y, X, Y, R, T] = CreateGrid(S.ImCubeUnProb{1}, 1./S.ppl0);
            bMaskSc = R >= RminSc & R <= RmaxSc;
            if ~isempty(TminSc) && ~isempty(TmaxSc),
                bMaskSc = bMaskSc & ...
                    ( (T >= TminSc & T <= TmaxSc) | (T >= mod2pi(TminSc + pi) & T <= mod2pi(TmaxSc + pi)) );
                
                %figure, imageschcit(x, y, bMaskSc), title('Check Scoring Region Mask')
            end

            % resample throughput data to pixels in scoring region
            if ~isempty(S.Sthpt),
                if S.runnum == 603, % SPC disc
                    error('SPC disc throughput data format needs to be reworked');
                    %                     dxyfudge = 0.9;
                    %                     rthpt = mean([Sth.src_dy(:).*dxyfudge./Sth.mmperlamD Sth.src_dx(:).*dxyfudge./Sth.mmperlamD],2);
                    %                     thpt = mean([Sth.thpt_dy(:) Sth.thpt_dx(:)],2);
                    %                     % normalize to clear region: lam/D > 8.5 & lam/D < 16.2
                    %                     thpt_norm = mean(thpt(abs(rthpt) >= 8.5 & abs(rthpt) < 16.2));
                    %                     thpt = thpt./ thpt_norm;
                    %                     % rthpt, thpt is now a lookup table to throughput for
                    %                     % cross-section
                    %                     % make a matrix
                    %                     Thpt = zeros(size(bMaskSc));
                    %                     Thpt(bMaskSc) = interp1(rthpt, thpt, R(bMaskSc));
                end
                Finterp = scatteredInterpolant(S.Sthpt.fovx(:), S.Sthpt.fovy(:), S.Sthpt.thpt(:));
                Thpt = ones(size(bMaskSc)); % avoid divide by zero
                Thpt(bMaskSc) = Finterp(X(bMaskSc), Y(bMaskSc));

            else
                Thpt = ones(size(bMaskSc));
                
            end
            
            % note: don't like nonzeros(...) because pixels within the dark
            % hole that happen to have zero measured intesnity do not
            % contribute to the mean.
            Contrast.total_lam = zeros(1,S.NofW);
            for iwv = 1:S.NofW,
                % NI total control region
                Contrast.total_lam(iwv) = mean(nonzeros( S.ImCube(:,:,S.imgindex(iwv)).*S.bMask ));
                % NI score region
                Contrast.score_lam(iwv) = mean(S.ImCubeUnProb{iwv}(bMaskSc)); 
                % Contrast score region
                Contrast.contr_lam(iwv) = mean(S.ImCubeUnProb{iwv}(bMaskSc)./Thpt(bMaskSc));
                
            end
            Contrast.mean = mean(Contrast.contr_lam);

            % mean coherent and incoherent contrast
            
            if isempty(S.CohInt),
                S.ReadReducedCube;
            end

            % scoring region
            Contrast.inco_lam = zeros(1,S.NofW);
            Contrast.co_lam = zeros(1,S.NofW);
            for iwv = 1:S.Nlamcorr,
                %bPampz = squeeze(S.bPampzero(iwv,:,:));
                %bMaskUse = ~bPampz & bMaskSc;
                
                bMaskUse = bMaskSc;
                % total incoherent
                Contrast.inco_lam(iwv) = mean(S.IncInt{iwv}(bMaskUse)./Thpt(bMaskUse));
                % estimated incoherent (those pixels where E was estimated)
                
                % Mix incoherent (those pixels where sigtbw < 0, Eest = 0)
                
                % coherent
                Contrast.co_lam(iwv)   = mean(S.CohInt{iwv}(bMaskUse)./Thpt(bMaskUse));
            end
            Contrast.inco_mean = mean(Contrast.inco_lam);
            Contrast.co_mean = mean(Contrast.co_lam);
            
        end % GetContrast

        function [Excel, sheet] = ContrastReportExcel(S, varargin)
            
            if ~isequal(computer, 'PCWIN64'),
                Excel = [];
                sheet = [];
                warning('Must Use PCWIN To Open Excel');
                return
            end

            Rmin = CheckOption('RminSc', S.RminSc, varargin{:});
            Rmax = CheckOption('RmaxSc', S.RmaxSc, varargin{:});
            
            [Excel, Workbook, Sheets, sheet] = Exlopen;

            
            Exlsetval(sheet, {'A2','A6'}, cellstr(num2str(S.NKTcenter(:)/S.NM,'%.1f')));
            Exlsetval(sheet, 'A7', 'Mean');
            
            for irr = 1:length(Rmin),
                col = num2column(irr+1);
                Ctmp = S.GetContrast('RminSc',Rmin(irr),'RmaxSc',Rmax(irr));
                Exlsetval(sheet, [col '1'], [num2str(Rmin(irr),'%.2f') ' - ' num2str(Rmax(irr),'%.2f')]);
                Exlsetval(sheet, {[col '2'],[col '6']}, Ctmp.contr_lam(:));
                Exlsetval(sheet, [col '7'], Ctmp.mean);
                
            end

            Exlsetval(sheet, 'A8', S.Sthpt.ThptCal_fn);
            
        end % ContrastReportExcel
        
        function S = ReadReducedCube(S)
            % reads the primary hdu of the reduced data fits
            
            if isempty(S.Reduced_fn),
                fprintf('no reduced data');
                return
            end
            
            RedData = fitsread(S.Reduced_fn); % primary hdu
            [nr, nc, nslices] = size(RedData);
            % each images slice is nr x nc
            
            % Primary: reducedcube
            %          for each wave: incoherent +              (0)
            %                         nppair probe amps +      (1,2,3,4)
            %                         coherent (r, i) +         (5,6)
            %                         starting efield +    (r,i)(6,7)
            %                         optimal reg - starting + (r,i)
            %                         each bscan - starting + (r,i)
            %                         pampzero (r,i)
            
            %nbscan = 5;
            %numperwv = 1 + S.Nppair + 2 + 2 + 2*(1 + nbscan) + 2;
            
            S.bPampzero = false(S.Nlamcorr, nr, nc);
            
            if isempty(S.bMask),  S.ReadMaskCube; end
            for iwl = 1:S.Nlamcorr,

                % incoherent image
                S.IncInt{iwl} = RedData(:,:,0*S.Nlamcorr+iwl);

                % probe fields
                for ipr = 1:S.Nppair,
                    S.ProbeAmp{iwl,ipr} = RedData(:,:,ipr*S.Nlamcorr+iwl);
                end

                % coherent field
                S.CohInt{iwl} = RedData(:,:,(1+S.Nppair)*S.Nlamcorr+iwl).^2 + RedData(:,:,(1+S.Nppair+1)*S.Nlamcorr+iwl).^2;

                % should be mask.* sqrt(photcorr(ii,iwl)).*
                S.E_t(iwl,:,:)	= (RedData(:,:,(1+S.Nppair)*S.Nlamcorr+iwl)+1i*RedData(:,:,(1+S.Nppair+1)*S.Nlamcorr+iwl));
                S.E_m(iwl,:,:) 	= (RedData(:,:,(3+S.Nppair)*S.Nlamcorr+iwl)+1i*RedData(:,:,(3+S.Nppair+1)*S.Nlamcorr+iwl));

                %                 S.dE_ro(iwl,:,:) = (RedData(:,:,8*S.Nlamcorr+iwl)+1i*RedData(:,:,9*S.Nlamcorr+iwl));	% regularization optimal
                %                 S.dE_r1(iwl,:,:) = (RedData(:,:,10*S.Nlamcorr+iwl)+1i*RedData(:,:,11*S.Nlamcorr+iwl));	% regularization 1
                %                 % need to read mask first

                % pampzero 1 + Nppair + 2*(1 + 1 + 1 + Nbscan + 1)-1
                % pampzero is logical, imag and real are the same logical
                % pampzero = true for pixels where: eest == 0 OR iinc < -()
                %    OR eestcond < eestcondlim
                %                 isl = S.Nppair + 2*(1 + 1 + 1 + S.Nbscan + 1+1)-1;
                %                 S.bPampzero(iwl,:,:) = logical(RedData(:,:,isl * S.Nlamcorr+iwl)); % + 1i*RedData(:,:,(isl + 1)*S.Nlamcorr+iwl);

                % mean coherent and incoherent contrast
                S.Contrast_inco(iwl) = mean(nonzeros(S.IncInt{iwl}.*S.mdMask));
                S.Contrast_co(iwl)   = mean(nonzeros(S.CohInt{iwl}.*S.mdMask));
                
                % if require all probes is true, pixels where abs(CohInt) == 0,
                % are pixels where probe amp < 0 for at least one probe,
                % image at this pixel is considered all incoherent. Let's
                % separate the incoherent into two components:
                bProbNeg = S.mdMask .* ~( abs(S.CohInt{iwl}) > 0 );
                S.IncIntEst{iwl} = (~bProbNeg).*S.IncInt{iwl};
                S.IncIntMix{iwl} =   bProbNeg .*S.IncInt{iwl}; % = unprobed Image at these pixels
                       
            end % for each wl
            
        end % ReadReducedCube
        
        function S = ReadMaskCube(S)
            
            % 1st slice is mask of control region (md)
            % next nlamcorr * 2 slices are real and imag occulter transmission profile

            MaskCubeData = fitsread(S.Reduced_fn,'image',1);
            
            % 1st slice is mask of control region, i.e. model.band[ilam0].md
            S.mdMask = squeeze(MaskCubeData(:,:,1));
            S.bMask  = S.mdMask > 0;
            
            
        end % ReadMaskCube
        
        function S = ReadProbeCube(S)
                        
            % for each probe and wavelength (from efc.tbif.task.step, search
            % probecube)
            %                   1 ..   Nppair*Nlamcorr = model probe amplitude for each pair
            %   Nppair*Nlamcorr+1 .. 2*Nppair*Nlamcorr = model probe phase for each pair
            % 2*Nppair*Nlamcorr+1 .. 3*Nppair*Nlamcorr = measured (Iplus + Iminus) - I0            
            % 3*Npparr*Nlamcorr+1 .. 4*Nppair*Nlamcorr = measured Iplus - Iminus
            % 4*Nppair*Nlamcorr+1 .. 4*Nppair*Nlamcorr + Nlamcorr  = probe estimate residual
            % 4*Nppair*Nlamcorr+Nlamcorr+1 .. 4*Nppair*Nlamcorr+2*Nlamcorr = probe estimate condition #
            
            %finfo = fitsinfo(S.Reduced_fn);
            %ProbeKwds = finfo.Image(2).Keywords;
            ProbeData = fitsread(S.Reduced_fn,'image',2);

            for iwl = 1:S.Nlamcorr,
                for ip = 1:S.Nppair,
                    
                    S.ProbeModel{iwl, ip} = squeeze(...
                        ProbeData(:,:,0*S.Nlamcorr*S.Nppair+(ip-1)*S.Nlamcorr+iwl) ...
                        .* exp(1i*...
                        ProbeData(:,:,1*S.Nlamcorr*S.Nppair+(ip-1)*S.Nlamcorr+iwl)) ...
                        );
                    
                    S.ProbeMeasAmp{iwl, ip} = sqrt(squeeze( ProbeData(:,:,2*S.Nlamcorr*S.Nppair + (ip-1)*S.Nlamcorr+ iwl) ));
                    S.ProbeMeasCross{iwl, ip} = squeeze( ProbeData(:,:,3*S.Nlamcorr*S.Nppair + (ip-1)*S.Nlamcorr + iwl) );
                    
                end % for each probe pair                
                S.ProbeRes{iwl} = squeeze( ProbeData(:,:,4*S.Nlamcorr*S.Nppair + iwl) );
                S.ProbeCon{iwl} = squeeze( ProbeData(:,:,4*S.Nlamcorr*S.Nppair + S.Nlamcorr + iwl) );
            end % for each wave
            
        end % ReadProbeCube
        
        function S = ReadImageCube(S)
            
            % unprobed images for total contrast:
            % see hcim/ly/lyotserver.py,
            % wavelength is outer loop, dmsettings (probes) inner loop
            %total       = ddtotal(:,:,1:7:(NofW*7));
            % # of images in ImCube = NofW * (2*nppair + 1)
            
            S.ImCube = fitsread(S.Rundir_fn);
            
            % unprobed images, each wavelength:
            for iwv = 1:S.NofW,
                imnum = (iwv-1)*(2*S.Nppair+1) + 1;
                S.ImCubeUnProb{iwv} = squeeze(S.ImCube(:,:,imnum));
                for ipr = 1:S.Nppair,
                    imnum = (iwv-1)*(2*S.Nppair+1) + 1 + (2*ipr-1);
                    ImPrPlus = squeeze(S.ImCube(:,:,imnum));
                    ImPrMinus = squeeze(S.ImCube(:,:,imnum+1));
                    S.ImCubeDelProb{iwv,ipr} = ImPrPlus - ImPrMinus;
                    S.ImCubeSigProb{iwv,ipr} = 0.5*(ImPrPlus + ImPrMinus) - S.ImCubeUnProb{iwv};
                end
            end
            
            % Useful Derived Parameters                       
            
        end % ReadImageCube
        
        function S = ReadDMvCube(S)
            
            finfo = fitsinfo(S.Rundir_fn);
            % Primary hdu is unprobed images
            % one Image hdu per DM
            S.Ndm = length(finfo.Image);
            
            for dmnum = 1:S.Ndm,
                S.DMvCube{dmnum} = fitsread(S.Rundir_fn,'image',dmnum);
            end
            
        end % ReadDMv

        function [hfig, him, Im] = DisplayImCubeImage(S, imnum)
            if nargin < 2,
                disp('usage: S.DisplayImCubeImage(imnum)');
                error('Not enough input arguments');
            end
            
            if isempty(S.ImCube),
                S.ReadImageCube;
            end
            
            Im = abs(squeeze(S.ImCube(:,:,imnum)));
            
            hfig = figure;
            him = imageschcit(Im); axis image, colorbar
            title(['Iteration #' num2str(S.iter) ', Image #' num2str(imnum)])
            
        end % DisplayImCubeImage
        
        function [hfig, ha] = DisplayImCubeUnProb(S, varargin)
            
            if isempty(S.ImCube),
                S.ReadImageCube;
            end
            
            % default options and set requested options
            %  val = CheckOption(sOpt, valDefault, varargin)
            bPlotLog = CheckOption('bLog', false, varargin{:});
            dispXYlim = CheckOption('xylim', S.XYlimDefault, varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            climopt = CheckOption('clim', [], varargin{:});
            ilam = CheckOption('ilam', 1:S.NofW, varargin{:});
            haxuse = CheckOption('hax', [], varargin{:});
            
            [x, y, X, Y, R] = CreateGrid(S.ImCubeUnProb{1}, 1./S.ppl0);
            xlim = dispXYlim*[-1 1]; ylim = xlim;
            
            Nlam = length(ilam);
            
            if isempty(haxuse), hfig = figure_mxn(1,Nlam); else hfig = gcf; end
            for iwv = 1:Nlam,
                iwvpl = ilam(iwv);
                if isempty(haxuse),
                    ha(iwv) = subplot(1,Nlam,iwv);
                else
                    ha(iwv) = haxuse(iwv);
                    axes(haxuse(iwv));
                end
                
                if bPlotLog,
                    imageschcit(x, y, log10(abs(S.ImCubeUnProb{iwvpl}))), axis image,
                    colorbartitle('log_{10} Norm Intensity')
                else
                    imageschcit(x, y, S.ImCubeUnProb{iwvpl}), axis image,
                    colorbar
                end
               clim(iwv,:) = get(gca,'clim');
               
               set(gca,'xlim',xlim,'ylim',ylim)
               xlabel('\lambda/D'), ylabel('\lambda/D')
               titlestr = ['Iter #' num2str(S.iter)];
               if ~isempty(S.NKTcenter),
                   titlestr = [titlestr ', Wave ' num2str(S.NKTcenter(iwvpl)/S.NM) 'nm'];
               end
               title(titlestr)
                               
            end % for each wavelength and subplot

            % overlay circles if requested
            DrawCircles(ha, drawRadii);
            DrawThetaLines(ha, drawTheta, drawRadii);

            % set each image plot to the same clim
            % auto-clim, unless specific clim requested
            if isempty(climopt),
                cmaxsort = sort(clim(:),'descend');
                %set(ha,'clim', [min(clim(:)) cmaxsort(2)]);
                set(ha,'clim', [max(clim(:,1)) min(clim(:,2))])
            else
                set(ha,'clim',climopt);
            end

            % radial plot
            Nr = 128;
            re = linspace(0, dispXYlim, Nr+1)';
            S.IntRad = cell(S.NofW,1);
            legstr = cell(1,S.NofW);
            for iwv = 1:S.NofW,
                Itmp = zeros(Nr,1);
                for ir = 1:Nr,
                    Itmp(ir) = mean(S.ImCubeUnProb{iwv}(R > re(ir) & R <= re(ir+1)));
                end % for ir
                S.IntRad{iwv} = Itmp;
                legstr{iwv} = ['Wave #' num2str(iwv)];
            end % for iwv
        
            S.rplot = mean([re(1:end-1) re(2:end)],2); % radii midway between edges
            figure, semilogy(S.rplot, [S.IntRad{:}]), grid
            if ~isempty(drawRadii),
                ylim = get(gca,'ylim');
                hold on
                for irad = 1:length(drawRadii),
                    plot(drawRadii(irad)*[1 1], ylim, '--r')
                end
                hold off
            end

            xlabel('Radius (\lambda/D)')
            ylabel('Average Normalized Intensity')
            legend(legstr{:},'location','north')
            title(['Iter #' num2str(S.iter)])
            
        end % DisplayImCubeUnProb(S)
        
        function [hfig, ha] = DisplayImCubeSigProb(S, varargin)
            
            if isempty(S.ImCube),
                S.ReadImageCube;
            end
            
            % default options and set requested options
            %  val = CheckOption(sOpt, valDefault, varargin)
            bPlotLog = CheckOption('bLog', false, varargin{:});
            dispXYlim = CheckOption('xylim', S.XYlimDefault, varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            climopt = CheckOption('clim', [], varargin{:});
            ilam = CheckOption('ilam', 1:S.NofW, varargin{:});
            ipro = CheckOption('iprobe', 1:S.Nppair, varargin{:});
            haxuse = CheckOption('hax', [], varargin{:});
            
            [x, y, X, Y, R] = CreateGrid(S.ImCubeSigProb{1,1}, 1./S.ppl0);
            xlim = dispXYlim*[-1 1]; ylim = xlim;
            
            Nlam = length(ilam);
            Npro = length(ipro);
            
            if isempty(haxuse), 
                hfig = figure_mxn(Npro,Nlam);
            else
                if ~isequal(size(haxuse),[Npro, Nlam]),
                    error('input haxuse must be size Npro, Nlam');
                end
                hfig = gcf;
            end
            
            for iwv = 1:Nlam,
                for ipr = 1:Npro,
                    iwvpl = ilam(iwv);
                    iprpl = ipro(ipr);
                    
                    if isempty(haxuse),
                        ha(ipr,iwv) = subplot(Npro,Nlam,(ipr-1)*Nlam + iwv);
                    else
                        ha(ipr,iwv) = haxuse(ipr,iwv);
                        axes(haxuse(ipr,iwv));
                    end
                    
                    if bPlotLog,
                        error('log plot not implemented');
                        %imageschcit(x, y, log10(abs(S.ImCubeUnProb{iwvpl}))), axis image,
                        %colorbartitle('log_{10} Norm Intensity')
                    else
                        imageschcit(x, y, S.ImCubeSigProb{iwvpl, iprpl}), axis image,
                        colorbar
                    end
                    clim(ipr,iwv,:) = get(gca,'clim');
               
                    set(gca,'xlim',xlim,'ylim',ylim)
                    xlabel('\lambda/D'), ylabel('\lambda/D')
                    titlestr = ['Iter #' num2str(S.iter)];
                    if ~isempty(S.NKTcenter),
                        titlestr = [titlestr ', Wave ' num2str(S.NKTcenter(iwvpl)/S.NM) 'nm'];
                    end
                    titlestr = [titlestr ', Pr #' num2str(iprpl)];
                    title(titlestr)
                
                end % for each probe and subplot
                
            end % for each wavelength 
            
            % overlay circles on all the axes if requested
            DrawCircles(ha, drawRadii);
            DrawThetaLines(ha, drawTheta, drawRadii);

            %             % set each image plot to the same clim
            %             % auto-clim, unless specific clim requested
            %             if isempty(climopt),
            %                 cmaxsort = sort(clim(:),'descend');
            %                 %set(ha,'clim', [min(clim(:)) cmaxsort(2)]);
            %                 set(ha,'clim', [max(clim(:,1)) min(clim(:,2))])
            %             else
            %                 set(ha,'clim',climopt);
            %             end
            if ~isempty(climopt),
                set(ha,'clim',climopt);
            end
            
        end % DisplayImCubeSigProb
        
        function [hfig, hax] = DisplayIncInt(S, varargin)
            % [hfig, hax] = DisplayIncInt(S, varargin)
            % options:
            %    'drawradii'
            %    'bLog' (default = true)
            %    'clim'
            %    'hax'
            %    'type' (default = 'normal'), 'Est', 'Mix'
            
            if isempty(S.IncInt),
                S.ReadReducedCube;
            end

            
            % options:
            % %  val = CheckOption(sOpt, valDefault, varargin)
            bLog = CheckOption('bLog', true, varargin{:});
            dispXYlim = CheckOption('xylim', S.XYlimDefault, varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            clim = CheckOption('clim', [], varargin{:});
            haxuse = CheckOption('hax', [], varargin{:}); % put image on this axes
            IncIntType = CheckOption('type', 'normal', varargin{:});
            
            %%%% end options

            switch lower(IncIntType),
                case 'normal',
                    plIncInt = S.IncInt;
                case 'est',
                    plIncInt = S.IncIntEst;
                case 'mix',
                    plIncInt = S.IncIntMix;
                otherwise,
                    error(['invalid type: ' IncIntType]);
            end
            
            if bLog,
                pFun  = @(a) real(log10(a)) .* (-1).^imag(log(a)/pi);
                pClim = @(a) [-10 0.99*max(real(log10(a(:))))];
                cbartitle = 'log_{10} Norm Intensity';
            else
                pFun  = @(a) (a<0);
                pClim = @(a) AutoClim(a);
                cbartitle = 'Norm Intensity';
            end
            
            if isempty(haxuse), hfig = figure_mxn(1,S.Nlamcorr); else hfig = gcf; end
            
            xlim = dispXYlim*[-1 1]; ylim = xlim;

            if isempty(clim),
                clim = pClim(pFun([plIncInt{:}]));
            end
            
            [x, y] = CreateGrid(plIncInt{1}, 1./S.ppl0);
            % auto-scale
            %Agg = [
            for iwv = 1:S.Nlamcorr,

                if isempty(haxuse),
                    hax(iwv) = subplot(1, S.Nlamcorr, iwv);
                    colormap(gray)
                else
                    hax(iwv) = haxuse(iwv);
                    axes(haxuse(iwv));
                end
                
                him(iwv) = imageschcit(x, y, pFun(plIncInt{iwv})); axis image
                xlabel('\lambda / D')
                ylabel('\lambda / D')
                %caxis(clim);

                colorbartitle(cbartitle)                
                if ~isempty(S.NKTcenter), strlam = ['\lambda = ' num2str(S.NKTcenter(iwv)/S.NM) 'nm']; else strlam = ['Wave #' num2str(iwv)]; end
                title(['it#' num2str(S.iter) ', ' strlam])
            
            end
            set(hax,'xlim',xlim,'ylim',ylim,'clim',clim);
            
            DrawCircles(hax, drawRadii);
            
            %             % radial plot
            %             for iwv = 1:S.Nlamcorr,
            %                 [fovrplot, IncIntrad] = RadialMean(x, y, plIncInt{iwv}, 128);
            %                 % remove negative inc int values so we can use log plot
            %                 IncIntrad(IncIntrad < 1e-11) = 1e-11;
            %                 qplot{1,iwv} = fovrplot;
            %                 qplot{2,iwv} = IncIntrad;
            %                 legstr{iwv} = ' ';
            %                 if ~isempty(S.NKTcenter), legstr{iwv} = [num2str(S.NKTcenter(iwv)/S.NM) 'nm']; end
            %             end
            %             figure, semilogy(qplot{:});
            %             hax(end+1) = gca;
            %             grid on
            %             set(gca,'xlim',[0 1].*xlim)
            %             %set(gca,'ylim',clim)
            %
            %             xlabel('Radius (\lambda/D)')
            %             ylabel('Unmodulated (Norm. Int.)')
            %             title(['Inc Int, it#' num2str(S.iter)])
            %
            %             if ~isempty(drawRadii),
            %                 hold on
            %                 for irad = 1:length(drawRadii),
            %                     plot(drawRadii(irad)*[1 1], get(gca,'ylim'), '--r')
            %                 end
            %                 hold off
            %             end
            %
            %             % legend must go after everything is plotted
            %             legend(legstr{:},'Location','North')

        end % DisplayIncInt
                
        function [hfig, hax] = DisplayProbeAmp(S, varargin)
            if isempty(S.ProbeAmp),
                S.ReadReducedCube;
            end
            
            % check options
            clim = CheckOption('clim', [], varargin{:});
            dispXYlim = CheckOption('xylim', S.XYlimDefault, varargin{:});
            
            xlimlamD = dispXYlim;
            [x, y] = CreateGrid(S.ProbeAmp{1,1}, 1./S.ppl0);
            
            % ProbeAmp is from reduced data cube, Primary HDU
            % arrange subplots
            Nwide = max([S.Nppair, S.Nlamcorr]);
            Nhigh = min([S.Nppair, S.Nlamcorr]);
            hfig(1) = figure_mxn(Nhigh, Nwide);
            for ip = 1:S.Nppair,
                for iw = 1:S.Nlamcorr,
                    if Nwide == S.Nlamcorr,
                        ipl = (ip-1)*S.Nlamcorr+iw;
                    else
                        ipl = (iw-1)*S.Nppair+ip;
                    end
                    hax(ipl) = subplot(Nhigh, Nwide, ipl);
                    imageschcit(x, y, S.ProbeAmp{iw,ip}), axis image,
                    xlabel('Cam X (\lambda/D)'), ylabel('Cam Y (\lambda/D)')
                    colorbartitle('Amplitude (Normalized)'),
                    set(gca,'xlim',xlimlamD*[-1 1],'ylim',xlimlamD*[-1 1])
                    title(['Iteration # ' num2str(S.iter) ' Wave #' num2str(iw) ' Probe # ' num2str(ip)])
                    meanProbeh = mean( abs(S.ProbeAmp{iw,ip}(S.bMask)).^2 );
                    ht = text(-xlimlamD, xlimlamD, [' Mean Probe Int = ' num2str(meanProbeh, '%.1e')]);
                    set(ht,'VerticalAlignment','top')
                    set(ht,'Color',[1 1 1])
                    
                    axclim(ipl,:) = get(hax(ipl),'clim');
                    
                end
            end
            
            % adjust clim
            if isempty(clim),
                clim = [min(axclim(:)) max(axclim(:))];
            end
            set(hax, 'clim', clim);

            if ~any(S.bPampzero(:)), return, end
            
            % again, but only show points where bPampzero = false
            % i.e. pixels used by EFC
            hfig(2) = figure_mxn(Nhigh, Nwide);
            for ip = 1:S.Nppair,
                for iw = 1:S.Nlamcorr,
                    if Nwide == S.Nlamcorr,
                        ipl = (ip-1)*S.Nlamcorr+iw;
                    else
                        ipl = (iw-1)*S.Nppair+ip;
                    end
                    subplot(Nhigh, Nwide, ipl),
                    bPamp = squeeze(S.bPampzero(iw,:,:));
                    probeamptmp = S.ProbeAmp{iw,ip};
                    probeamptmp(bPamp) = 0.0;
                    imageschcit(x, y, probeamptmp), axis image, 
                    xlabel('Cam X (\lambda/D)'), ylabel('Cam Y (\lambda/D)')
                    colorbartitle('Amplitude (Normalized)'),
                    set(gca,'xlim',xlimlamD*[-1 1],'ylim',xlimlamD*[-1 1])
                    title(['Iteration # ' num2str(S.iter) ' Wave #' num2str(iw) ' Probe # ' num2str(ip) ' Amplitude'])
                    
                end
            end

        end % DisplayProbeAmp
        
        function [pl, hfig, hax] = DisplayOneProbeAmp(S, ip, varargin)
            if isempty(S.ProbeAmp),
                S.ReadReducedCube;
            end

            dispXYlim = CheckOption('xylim', S.XYlimDefault, varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
                
            % display one probe at all wavelengths and plot cross sections
            [x, y] = CreateGrid(S.ProbeAmp{1,1}, 1./S.ppl0);
            xlimlamD = 22*[-1 1];
            hfig = figure_mxn(1,S.Nlamcorr);
            for iw = 1:S.Nlamcorr,
                hax(iw) = subplot(1,S.Nlamcorr,iw);
                imageschcit(x, y, (S.ProbeAmp{iw, ip}).^2), axis image, colorbartitle('Norm Intensity')
                set(gca,'xlim',xlimlamD,'ylim',xlimlamD)
                title(['Iteration # ' num2str(S.iter) ' Wave #' num2str(iw) ' Probe # ' num2str(ip) ' Amplitude'])
                                
                pl{1,iw} = [x(:) y(:)];
                pl{2,iw} = [S.ProbeAmp{iw, ip}(y==0, :).' S.ProbeAmp{iw, ip}(:, x==0)].^2;
                
            end 
            DrawCircles(hax, drawRadii);
            
            figure, semilogy(pl{:}), grid
            set(gca,'xlim',xlimlamD)
            xlabel('X, Y (\lambda/D)')
            ylabel('Norm Intensity')
            title(['Iter # ' num2str(S.iter) ', Probe # ' num2str(ip) ])
            
            
        end % DisplayOneProbeAmp
        
        function [hfig, ha] = DisplayProbeCube(S, iwvplot, varargin)
            % hfig = DisplayProbeCube(S, iwvplot)
            % ProbeCube is 2nd HDU of reduced data cube
            if isempty(S.ProbeModel),
                S.ReadProbeCube;
            end
            
            if ~exist('iwvplot','var') || isempty(iwvplot),
                warning('wavelength # not specified, plotting wave #1');
                iwvplot = 1;
            end
            
            dispXYlim = CheckOption('xylim', S.XYlimDefault, varargin{:});
            bLog = CheckOption('blog', true, varargin{:});

            [x, y] = CreateGrid(S.ProbeModel{1,1}, 1./S.ppl0);
            
            ip = 1;
            figure, plot(y, [abs(S.ProbeModel{iwvplot,ip}(:,x==0)), ...
                    S.ProbeMeasAmp{iwvplot,ip}(:,x==0)].^2), grid
            xlabel('Y (\lambda / D)')
            legend(['it ' num2str(S.iter) ' model'],['it ' num2str(S.iter) ' Measured'])
            title(['Cross Section, it ' num2str(S.iter) ', Wave #' num2str(iwvplot) ', Measure Probe #' num2str(ip)])

            if bLog,
                ImModel = real(log10(abs(S.ProbeModel{iwvplot,ip}).^2));
                ImMeasure = real(log10(S.ProbeMeasAmp{iwvplot,ip}.^2));
                sctitle = 'log_{10} Norm Intensity';
            else
                ImModel = abs(S.ProbeModel{iwvplot,ip}).^2;
                ImMeasure = S.ProbeMeasAmp{iwvplot,ip}.^2;
                sctitle = 'Linear Norm Intensity';
            end
                
            hfig = figure_mxn(2,S.Nppair);
            for ip = 1:S.Nppair,
                ha(ip) = subplot(2,S.Nppair,ip);
                imageschcit(x,y,ImModel), axis image, colorbartitle(sctitle)
                set(gca,'xlim',dispXYlim*[-1 1],'ylim',dispXYlim*[-1 1])
                xlabel('\lambda / D'), ylabel('\lambda / D')
                title(['wave #' num2str(iwvplot) ', it ' num2str(S.iter) ', Model Probe #' num2str(ip)])
                
                ha(S.Nppair+ip) = subplot(2,S.Nppair,S.Nppair+ip);
                imageschcit(x,y,ImMeasure), axis image, colorbartitle(sctitle)
                set(gca,'xlim',dispXYlim*[-1 1],'ylim',dispXYlim*[-1 1])
                xlabel('\lambda / D'), ylabel('\lambda / D')
                title(['wave #' num2str(iwvplot) ', it ' num2str(S.iter) ', Measure Probe #' num2str(ip)])
            end

        end % DisplayProbeCube
        
        function [hfig, hax] = DisplayCohInt(S, varargin)
            % [hfig, hax] = DisplayCohInt(S, varargin)
            % options:
            %    'drawradii'
            %    'bLog' (default = true)
            %    'clim'
            %    'hax'
            
            if isempty(S.E_t),
                S.ReadReducedCube;
            end
            
            for iwv = 1:S.Nlamcorr,
                CohInt{iwv} = abs(squeeze(S.E_t(iwv,:,:))).^2;
            end

            % options:
            % %  val = CheckOption(sOpt, valDefault, varargin)
            bLog = CheckOption('blog', true, varargin{:});
            dispXYlim = CheckOption('xylim', S.XYlimDefault, varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            clim = CheckOption('clim', [], varargin{:});
            haxuse = CheckOption('hax', [], varargin{:}); % put image on this axes
            
            %%%% end options
            
            if bLog,
                pFun  = @(a) real(log10(a));
                pClim = @(a) [-10 0.99*max(real(log10(a(:))))];
                cbartitle = 'log_{10} Norm Intensity';
            else
                pFun  = @(a) (a);
                pClim = @(a) AutoClim(a);
                cbartitle = 'Norm Intensity';
            end
            
            if isempty(haxuse), hfig = figure_mxn(1,S.Nlamcorr); else hfig = gcf; end
            
            xlim = dispXYlim*[-1 1]; ylim = xlim;

            if isempty(clim),
                clim = pClim(pFun([CohInt{:}]));
            end
            
            [x, y] = CreateGrid(CohInt{1}, 1./S.ppl0);

            for iwv = 1:S.Nlamcorr,
                if isempty(haxuse),
                    hax(iwv) = subplot(1, S.Nlamcorr, iwv);
                else
                    hax(iwv) =  haxuse(iwv);
                    axes(haxuse(iwv));
                end
              
                him(iwv) = imageschcit(x, y, pFun(CohInt{iwv})); axis image
                xlabel('\lambda / D')
                ylabel('\lambda / D')
                %caxis(clim);
                colorbartitle(cbartitle)                
                if ~isempty(S.NKTcenter), strlam = ['\lambda = ' num2str(S.NKTcenter(iwv)/S.NM) 'nm']; else strlam = ['Wave #' num2str(iwv)]; end
                title(['it#' num2str(S.iter) ', ' strlam])
            
            end
            set(hax,'xlim',xlim,'ylim',ylim,'clim',clim);
            
            DrawCircles(hax, drawRadii);
            %
            %             % radial plot
            %             for iwv = 1:S.Nlamcorr,
            %                 [fovrplot, CohIntrad] = RadialMean(x, y, S.CohInt{iwv}, 128);
            %                 % remove negative inc int values so we can use log plot
            %                 IncIntrad(CohIntrad < 1e-11) = 1e-11;
            %                 qplot{1,iwv} = fovrplot;
            %                 qplot{2,iwv} = CohIntrad;
            %                 legstr{iwv} = ' ';
            %                 if ~isempty(S.NKTcenter), legstr{iwv} = [num2str(S.NKTcenter(iwv)/S.NM) 'nm']; end
            %             end
            %             figure, semilogy(qplot{:});
            %             hax(end+1) = gca;
            %             grid on
            %             set(gca,'xlim',[0 1].*xlim)
            %             %set(gca,'ylim',clim)
            %
            %             xlabel('Radius (\lambda/D)')
            %             ylabel('Modulated (Norm. Int.)')
            %             title(['Mod Int, it#' num2str(S.iter)])
            %
            %             if ~isempty(drawRadii),
            %                 hold on
            %                 for irad = 1:length(drawRadii),
            %                     plot(drawRadii(irad)*[1 1], get(gca,'ylim'), '--r')
            %                 end
            %                 hold off
            %             end
            %
            %             % legend must go after everything is plotted
            %             legend(legstr{:},'Location','North')
            
            
        end % DisplayCohInt
        
        function [hfig, haxlist] = DisplayAllInt(S, varargin)
            % display large table of unprobed, coh int, inc int images
            
            % defaults that might be different
            varargin{end+1} = 'bLog'; varargin{end+1} = true;
            
            hfig = figure_mxn(3,S.Nlamcorr);
            haxlist = zeros(3,S.Nlamcorr);
            
            % unprobed images
            figure(hfig);
            for ii = 1:S.Nlamcorr,
                haxlist(1,ii) = subplot(3,S.Nlamcorr,ii);
            end
            S.DisplayImCubeUnProb('hax',haxlist(1,:),varargin{:});

            % Coh Int
            figure(hfig);            
            for ii = 1:S.Nlamcorr,
                haxlist(2,ii) = subplot(3,S.Nlamcorr,ii+S.Nlamcorr);
            end

            S.DisplayCohInt('hax',haxlist(2,:), varargin{:});

            % Inc Int
            figure(hfig);            
            for ii = 1:S.Nlamcorr,
                haxlist(3,ii) = subplot(3,S.Nlamcorr,ii+2*S.Nlamcorr);
            end
            S.DisplayIncInt('hax',haxlist(3,:), varargin{:});

            % make common clim
            %climlist = get(haxlist,'clim');
            %set(haxlist,'clim',[min([climlist{:}]) max([climlist{:}])])

            set(haxlist,'clim',[-9 -6.5])

        end
        
        function [hfig, haxlist] = DisplayIncCohInt(S, varargin)
            %    create large display of:
            %      CohInt
            %      IncInt Estimated
            %      IncInt Mix (not estimated)
            %
            %    options:
            %      'hfig', hfig, display to hfig figure window
            
            % defaults that might be different
            varargin{end+1} = 'bLog'; varargin{end+1} = true;
            
            hfig = CheckOption('hfig', [], varargin{:});
            
            if isempty(hfig),
                hfig = figure_mxn(3,S.Nlamcorr);
            else
                figure_mxn(hfig, 3, S.Nlamcorr);
            end
            
            haxlist = zeros(3,S.Nlamcorr);
            
            % estimated coherent images
            figure(hfig);
            for ii = 1:S.Nlamcorr,
                haxlist(1,ii) = subplot(3,S.Nlamcorr,ii);
            end
            S.DisplayCohInt('hax',haxlist(1,:),varargin{:});

            % Incoherent estimated
            figure(hfig);            
            for ii = 1:S.Nlamcorr,
                haxlist(2,ii) = subplot(3,S.Nlamcorr,ii+S.Nlamcorr);
            end

            S.DisplayIncInt('type','est','hax',haxlist(2,:), varargin{:});

            % Inc Int Mix not estimated
            figure(hfig);            
            for ii = 1:S.Nlamcorr,
                haxlist(3,ii) = subplot(3,S.Nlamcorr,ii+2*S.Nlamcorr);
            end
            S.DisplayIncInt('type','mix','hax',haxlist(3,:), varargin{:});

            % make common clim
            %climlist = get(haxlist,'clim');
            %set(haxlist,'clim',[min([climlist{:}]) max([climlist{:}])])

            set(haxlist,'clim',[-9 -6.5])

            % Now Add Labels to Each Row!!!
            posrow1 = get(haxlist(1,1),'Position');
            haxlab1 = axes('Position',[0.07 posrow1(2)+0.5*posrow1(4) 0.01 0.01]);
            axis off
            htlab(1) = text(haxlab1,0,0,'Estimated Modulated');
            
            posrow2 = get(haxlist(2,1),'Position');
            haxlab2 = axes('Position',[0.07 posrow2(2)+0.5*posrow2(4) 0.01 0.01]);
            axis off
            htlab(2) = text(haxlab2,0,0,'Estimated Unmodulated');
            
            posrow3 = get(haxlist(3,1),'Position');
            haxlab3 = axes('Position',[0.07 posrow3(2)+0.5*posrow3(4) 0.01 0.01]);
            axis off
            htlab(3) = text(haxlab3,0,0,'|Probe Amp| < 0 (no estimate)');
            
            set(htlab ...
                ,'HorizontalAlignment','center' ...
                ,'Rotation',90 ...
                ,'FontSize',14 ...
                ,'Color','r' ...
                );
%             
        end % DisplayIncCohInt
        
        function [hfig, ha] = DisplayEfields(S, iwvplot)
            if ~exist('iwvplot','var') || isempty(iwvplot),
                warning('wavelength # not specified, plotting wave #1');
                iwvplot = 1;
            end
            
            if isempty(S.E_t),
                S.ReadReducedCube;
            end
            
            xlimlamD = 22*[-1 1];
            sRI = ['run #' num2str(S.runnum) ', iter #' num2str(S.iter)];
            
            E_t = squeeze(S.E_t(iwvplot,:,:));
            E_m = squeeze(S.E_m(iwvplot,:,:));

            hfig = figure_mxn(2,3);
            [x, y] = CreateGrid(E_t, 1./S.ppl0);
            
            ha = [];
            ha(end+1) = subplot(2,3,1); imageschcit(x,y,abs(E_t)), axis image, colorbar,
            title([sRI ', |E_t|, wave #' num2str(iwvplot)])
            set(gca,'xlim',xlimlamD,'ylim',xlimlamD);
            ha(end+1) = subplot(2,3,2); imageschcit(x,y,real(E_t)), axis image, colorbar,
            title([sRI ', real(E_t), wave #' num2str(iwvplot)])
            set(gca,'xlim',xlimlamD,'ylim',xlimlamD);
            ha(end+1) = subplot(2,3,3); imageschcit(x,y,imag(E_t)), axis image, colorbar,
            title([sRI ', imag(E_t), wave #' num2str(iwvplot)])
            set(gca,'xlim',xlimlamD,'ylim',xlimlamD);
            
            % equalize the clim range
            E_t_b = E_t(S.bMask);
            set(ha(1),'clim', AutoClim(abs(E_t_b)))
            set(ha(2:3), 'clim', AutoClim([real(E_t_b); imag(E_t_b)],'two-sided','symmetric'))

            ha(end+1) = subplot(2,3,4); imageschcit(x,y,abs(E_m)), axis image, colorbar,
            title([sRI ', |E_m|, wave #' num2str(iwvplot)])
            set(gca,'xlim',xlimlamD,'ylim',xlimlamD);
            ha(end+1) = subplot(2,3,5); imageschcit(x,y,real(E_m)), axis image, colorbar,
            title([sRI ', real(E_m), wave #' num2str(iwvplot)])
            set(gca,'xlim',xlimlamD,'ylim',xlimlamD);
            ha(end+1) = subplot(2,3,6); imageschcit(x,y,imag(E_m)), axis image, colorbar,
            title([sRI ', imag(E_m), wave #' num2str(iwvplot)])
            set(gca,'xlim',xlimlamD,'ylim',xlimlamD);

            % equalize the clim range
            E_m_b = E_m(S.bMask);
            set(ha(4),'clim', AutoClim(abs(E_m_b)))
            set(ha(5:6), 'clim', AutoClim([real(E_m_b); imag(E_m_b)], 'two-sided','symmetric'))

            for iha = 1:length(ha), xlabel(ha(iha), '\lambda/D'), ylabel(ha(iha),'\lambda/D'), end
            
        end % DisplayEfields
        
        function [hfig, ha] = DisplayDEfields(S, Sref, hfig)
            
            if isempty(S.E_t),
                S.ReadReducedCube;
            end

            if isempty(Sref.E_t),
                Sref.ReadReducedCube;
            end

            if isempty(S.ImCube),
                S.ReadImageCube;
            end
            
            if isempty(Sref.ImCube),
                Sref.ReadImageCube;
            end
            
            [nw, nr, nc] = size(S.E_t);
            
            xlimlamD = 22*[-1 1];
            sRI = ['run #' num2str(S.runnum) ', iter #' num2str(S.iter) '--' num2str(Sref.iter)];
            
            %             E_t = squeeze(S.E_t(iwvplot,:,:));
            %             E_m = squeeze(S.E_m(iwvplot,:,:));

            % top row = real(DE_t)
            % middle row = imag(DE_t)
            % bottom row = UnProbInt - UnProbInt(ref)
            Nplr = 4;
            if exist('hfig','var') && isa(hfig,'matlab.ui.Figure'),
                figure(hfig)
            else,
                hfig = figure_mxn(Nplr,S.NofW);
            end
            [x, y] = CreateGrid([nc nr], 1./S.ppl0);
            ha = zeros(Nplr,S.NofW);
            %
            dE = S.E_t - Sref.E_t;
            
            % FFT of dE
            Nlam2 = 1200; % a guess right now, about ppl0*D in pixels
            Nfft = Nlam2 ./(mean(S.NKTcenter)./S.NKTcenter);
            for iwv = 1:S.NofW,
                dEtmp = squeeze(dE(iwv,:,:));
                % need to eliminate outliers
                rmsAmp = rms(dEtmp(abs(dEtmp)>0));
                dEtmp(abs(dEtmp)>3*rmsAmp) = 0;

                dEtmp = PadImArray( dEtmp, Nfft(iwv) );
                Epup{iwv} = fftshift(fft2(fftshift(dEtmp)));
            end
            
            climE = zeros(2*S.NofW,2);
            climI = zeros(S.NofW,2);
            for iwv = 1:S.NofW,
                ipu = iwv+0*S.NofW;
                ipr = iwv+1*S.NofW;
                ipi = iwv+2*S.NofW;
                ipf = iwv+3*S.NofW;
                
                % change in unprobed image intensity
                ha(1,iwv) = subplot(Nplr, S.NofW, ipu);
                him = imageschcit(x,y, S.ImCubeUnProb{iwv}-Sref.ImCubeUnProb{iwv}); colorbar
                title([sRI ', \Delta UnPr Int, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])
                climI(iwv,:) = AutoClim(get(him,'CData'),'symmetric');

                ha(2,iwv) = subplot(Nplr, S.NofW, ipr);
                him = imageschcit(x,y,squeeze(real(dE(iwv,:,:)))); colorbar
                title([sRI ', real{\DeltaE}, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])
                climE(ipr,:) = AutoClim(get(him,'CData'),'symmetric');
                
                ha(3,iwv) = subplot(Nplr, S.NofW, ipi);
                him = imageschcit(x,y,squeeze(imag(dE(iwv,:,:)))); colorbar
                title([sRI ', imag{\DeltaE}, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])
                climE(ipi,:) = AutoClim(get(him,'CData'),'symmetric');

                ha(4,iwv) = subplot(Nplr, S.NofW, ipf);
                [xf, yf] = CreateGrid(Epup{iwv});
                him = imageschcit(xf, yf, abs(Epup{iwv}));
                xfyflim = 4*46.3;
                set(gca,'xlim',xfyflim*[-1 1],'ylim',xfyflim*[-1 1])
                title([sRI ', FFT{\DeltaE}, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])
                
            end
            set(ha(1,:),'clim',median(climI))
            set(ha(2:3,:),'clim',median(climE))
            set(ha(1:3,:),'xlim',xlimlamD,'ylim',xlimlamD);
            
        end % DisplayDEfields

        function DisplayPampzero(S)
            
            [nw, ny, nx] = size(S.bPampzero);
            
            [x, y] = CreateGrid([nx ny], 1./S.ppl0);
            xlim = dispXYlim*[-1 1]; ylim = xlim;

            hfig = figure_mxn(1,S.NofW);
            for iw = 1:S.NofW,
               ha(iw) = subplot(1,S.NofW,iw);
               
               imageschcit(x, y, squeeze(S.bPampzero(iw,:,:))), axis image
               xlabel('\lambda / D')
               ylabel('\lambda / D')
               
               set(gca,'xlim',xlim,'ylim',ylim)
               
            end
            
        end % DisplayPampzero

        function [hfig, hax] = DisplayDMv(S, varargin)
            % [hfig, hax] = S.DisplayDMv
            % [hfig, hax] = S.DisplayDMv(Sref)
            % [hfig, hax] = S.DisplayDMv(refDM1v_fits, refDM2v_fits)
            
            if isempty(S.DMvCube)
                S.ReadDMvCube;
            end

            % extract the DV v to plot
            for idm = 1:S.Ndm,
                DMvtmp = squeeze(S.DMvCube{idm}(:,:,1));
                rmsDMv(idm) = rms(DMvtmp(DMvtmp>0));
                
                DMv{idm} = DMvtmp;
                strDM{idm} = ['it#' num2str(S.iter) ...
                    ', DM' num2str(idm) ...
                    ', ' num2str(rmsDMv(idm),'%.4f') 'V rms'];
            end

            switch length(varargin),
                case 0,
                    refDMv = [];
                    strRefDM = [];

                case 1,
                    if ~isa(varargin{1},'CRunData'),
                        error('usage');
                    end
                    
                    Sref = varargin{1};
                    if isempty(Sref.DMvCube)
                        Sref.ReadDMvCube;
                    end
                    
                    for idm = 1:Sref.Ndm,
                        refDMv{idm} = squeeze(Sref.DMvCube{idm}(:,:,1));
                        strRefDM{idm} = ['it#' num2str(Sref.iter) ', DM' num2str(idm)];
                    end
                    
                case 2,
                    [refDM1v_fn, refDM2v_fn] = deal(varargin{1:2});
                    
                    refDMv{1} = fitsread(PathTranslator(refDM1v_fn));
                    aatmp = regexp(refDM1v_fn, '/', 'split');
                    strRefDM{1} = pwd2titlestr(aatmp{end});
                                                    
                    refDMv{2} = fitsread(PathTranslator(refDM2v_fn));
                    aatmp = regexp(refDM2v_fn, '/', 'split');
                    strRefDM{2} = pwd2titlestr(aatmp{end});


                otherwise %else
                    error('usage');
                    
            end % switch nargin

            if isempty(refDMv),
                Nr = 1;
            else
                Nr = 2;
            end
            hfig = figure_mxn(Nr,S.Ndm);
            
            for idm = 1:S.Ndm,

                % plot S DMv
                hax(idm) = subplot(Nr, S.Ndm, idm);
                imageschcit(DMv{idm}), 
                colorbartitle('Vmu')
                title(strDM{idm})
                
                % if Ref defined:
                if ~isempty(refDMv),
                    hax(idm+S.Ndm) = subplot(Nr, S.Ndm, idm+S.Ndm);
                
                    dDMv = DMv{idm} - refDMv{idm};
                    rmsdDMv = rms(dDMv(abs(dDMv)>0));

                    imageschcit(dDMv)
                    colorbartitle('Vmu')
                    title(['\Delta ' strRefDM{idm} ', ' num2str(rmsdDMv,'%.4f') 'V rms'])
                end
                
            end % for idm
            
            %             fprintf('rms dDMv1 = %.3f Vmu\n',rmsdDMv1);
            %             fprintf('rms dDMv2 = %.3f Vmu\n',rmsdDMv2);
            %
        end % DisplayDMv
        
        %
        %         function DisplayRadialContrast(S)
        %
        %         end % DisplayRadialContrast
    
        function [hfig, hax] = DisplayDelProbes(S, iwv, varargin)
            
            if isempty(S.ImCube),
                S.ReadImageCube;
            end
            
            if ~exist('iwv','var') || isempty(iwv),
                iwv = 1;
                warning('wave # not specified, displaying wave #1');
            end
            
            hfig = figure_mxn(1,3);
            for ii = 1:3,
                hax(ii) = subplot(1,3,ii);
                imageschcit(real(log10(S.bMask.*S.ImCubeDelProb{iwv,ii}))), colorbar,
            end
        end % DisplayDelProbes
% 
%         % 
%         S.ProbeRes
    end % methods
        
end % classdef

% utilities
function DrawCircles(hax, drawRadii)

if isempty(drawRadii),
    return
end

hax = hax(:);

for iax = 1:length(hax),
    axes(hax(iax));
    hold on
    for irad = 1:length(drawRadii),
        draw_circle([0 0], 2*drawRadii(irad), 1, 'r');
    end
    hold off
end % for each axes

end % DrawCircles

function DrawThetaLines(hax, drawTheta, drawRadii)

if isempty(drawTheta),
    return
end

hax = hax(:);

for iax = 1:length(hax),
    axes(hax(iax));

    if isempty(drawRadii),
        r0 = 0;
        r1 = max(get(hax(iax),'xlim'));
    else
        r0 = min(drawRadii);
        r1 = max(drawRadii);
    end
    

    
    hold on
    for ith = 1:length(drawTheta),
        th = drawTheta(ith) + pi/2;
        plot( [r0 r1]*cos(th), [r0 r1]*sin(th), '-r');
        plot(-[r0 r1]*cos(th),-[r0 r1]*sin(th), '-r');
        
    end % for each theta
    hold off;

end % for each axes

end % DrawThetaLines
