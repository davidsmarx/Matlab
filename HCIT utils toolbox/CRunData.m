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
% [hfig, hax] = DisplayDMvProbe(S)
%
% [hfig, ha] = DisplayDEfields(S, Sref, hfig)
%

    properties
        
        % SPC disc:
        lambda          = [542 553 565 577 588]*CConstants.NM;
        %ilamcorr        = [0, 1, 2, 3, 4];
        ppl0            
        PIAAMAG         = 1;   % only for PIAA testbed

        RminSc          = 6.3; % system (back-end) lam/D
        RmaxSc          = 19.5;
        ThminSc         = [];
        ThmaxSc         = [];
        YminSc          = -Inf;
        YmaxSc          = Inf;
        XminSc          = -Inf;
        XmaxSc          = Inf;
        bscan           = [];
        Nbscan          = 0;
        betamin
        Results_pn = '';
        Rundir_pn  = 'rundir/';         % always relative to Results_pn
        Reduced_pn = 'rundir/reduced/'; % always relative to Results_pn
        PCtemp_pn  = 'C:\Users\dmarx\HCITdatatemp\';
        S383temp_pn= '';
        Rundir_fn  = ''; % filenames include full path assembled in the init 
        Reduced_fn = '';
        svd_fn = '';
        sSVD
        
        runnum
        iter      = 0;
        timestamp   % taken from raw camera image file date
        
        % results:
        IncInt
        IncIntEst   % pixels where inc int < 0, inc int is fixed to = eps
        IncIntMix   % part of UnProbed Image where any probeamp <= 0
        IncIntFullBand 
        IncIntEstFullBand % mean across all subbands

        CohInt
        CohIntEst   % pixels where inc int < 0, coh int is fixed to = unprobed
        CohIntFullBand % mean across all subbands
        E_t
        E_m
        dE_ro
        dE_r1
        bPampzero
        dE_optimal
        dE_bscan
        bMask_badpix

        ProbeAmp % measured
        
        % control region and score region masks
        mdMask
        bMask
        bMaskSc
        
        ImCube
        ImCubeUnProb    % {iwv}
        ImCubeUnProbFullBand % mean across all subbands
        ImCubeDelProb   % {iwv,ipr}
        ImCubeSigProb   % {iwv,ipr} = sigtbw  in tbif.task.probes
        ImCubeContrast  % {iwv} = ImCubUnProb / Thpt
        ImCubeCohContrast
        ImCubeIncContrast
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
        NormIntensity_total
        NormIntensity_inco
        NormIntensity_co
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
        XlimDefault = []; % lam/D
        YlimDefault = []; % lam/D
        DrawradiiDefault = [];
        DrawthetaDefault = [];

        Sppt;
        
        debug = false;
        
    end % properties
        
    % reduced fits:
    % Primary: reducedcube
    %   Image: maskcube
    %   Image: probecube
    
    methods
        function S = CRunData(runnum, iter, sOptin, varargin)
            % S = CRunData(runnum, iter, sOptin, varargin)
            %
            % sOptin = struct, overides default class member values
            % varargin = list of methods to execute

            if nargin == 0,
                % return an empty instance
                return
            end
            
            S.runnum = runnum;
            S.iter = iter;

            % default paths based on runnum
            % 'rundir_pn' and 'reduced_pn' are always relative to
            % Results_pn
            switch S.runnum,
                case 0 % DST
                    S.Results_pn = '/home/dmarx/ln_dst_data/hcim/EFC/HLC/run000/';
                    S.XlimDefault = 12 * [-1 1];
                    S.YlimDefault = 12 * [-1 1];
                    
                case 001,
                    S.Results_pn = '/home/dmarx/HCIT/MCB/hcim_model2_run001/results/run001/';
                    
                    S.XlimDefault = 10 * [-1 1];
                    S.YlimDefault = 10 * [-1 1];
                    S.DrawradiiDefault = [3.0 9.0];

                case 10 % DST with BMC50.B DM at dm1
                    S.Results_pn = '/home/dmarx/ln_dst_data/EFC/HLC/run010/';
                    S.S383temp_pn= '/home/dmarx/HCIT/DST/hcim_testbed_run010/results/';
                    S.XlimDefault = 12 * [-1 1];
                    S.YlimDefault = 12 * [-1 1];
                    S.DrawradiiDefault = [3.0 9.0];
                    S.DrawthetaDefault = 180*[-0.5 0.5]*CConstants.P;
                    
                    S.RminSc    = 3.0; % lam/D
                    S.RmaxSc    = 9.0;
                    S.ThminSc   = (90 - 90)*CConstants.P;
                    S.ThmaxSc   = (90 + 90)*CConstants.P;

                    % overwritten if camera image is found
                    S.NKTupper = [533.5, 555.5, 577.5]*S.NM;
                    S.NKTlower = [522.5, 544.5, 566.5]*S.NM;
                    S.NKTcenter = mean([S.NKTupper; S.NKTlower]);
                    
                case 12 % DST with BMC50.B DM at dm1 & BMC 50.A dm2
                    S.Results_pn = '/home/dmarx/ln_dst_data/EFC/HLC/run012/';
                    S.S383temp_pn= '/home/dmarx/HCIT/DST/hcim_testbed_run012/results/';
                    S.XlimDefault = 12* [-1 1];
                    S.YlimDefault = 12* [-1 1];
                    S.DrawradiiDefault = [3.0 9.0];
                    
                    S.RminSc    = 3.0; % lam/D
                    S.RmaxSc    = 9.0;

                    % overwritten if camera image is found
                    S.NKTupper = [533.5, 555.5, 577.5]*S.NM;
                    S.NKTlower = [522.5, 544.5, 566.5]*S.NM;
                    S.NKTcenter = mean([S.NKTupper; S.NKTlower]);

                    S.ppl0 = 4.52;
                    
                case 13 % DST model2
                    S.Results_pn = '/home/dmarx/HCIT/DST/hcim_model2_run013/results/run013/';
                    S.S383temp_pn = S.Results_pn; % it's all local
                    
                    S.XlimDefault = 12 * [-1 1];
                    S.YlimDefault = 12 * [-1 1];

                    S.DrawradiiDefault = [3.0 9.0];
                    
                    S.RminSc    = 3.0; % lam/D
                    S.RmaxSc    = 9.0;

                    % overwritten if camera image is found
                    S.NKTupper = [533.5, 555.5, 577.5]*S.NM;
                    S.NKTlower = [522.5, 544.5, 566.5]*S.NM;
                    S.NKTcenter = mean([S.NKTupper; S.NKTlower]);

                    S.ppl0 = 4.45;

                case 100 % PIAA
                    S.Results_pn = '/proj/piaacmc/EFC/data/run100/';
                    S.S383temp_pn= '/home/dmarx/HCIT/PIAA/hcim_testbed_run100/results/';
                    
                    S.XlimDefault = 12 * [-1 1];
                    S.YlimDefault = 12 * [-1 1];

                    S.PIAAMAG = 1.12; % should get this from config
                    S.DrawradiiDefault = S.PIAAMAG*[1.8 9.0];
                    
                    % config_piaa_20210524.py
                    S.RminSc    = S.PIAAMAG * 4.0; % back-end (system) lam/D
                    S.RmaxSc    = S.PIAAMAG * 9.0;
                    S.YmaxSc    = S.PIAAMAG *-1.8;
                    S.YminSc    = -Inf;

                    % overwritten if camera image is found
                    S.NKTupper = [628.6]*S.NM; %[533.5, 555.5, 577.5]*S.NM;
                    S.NKTlower = [641.3]*S.NM; %[522.5, 544.5, 566.5]*S.NM;
                    S.NKTcenter = mean([S.NKTupper; S.NKTlower]);

                    %S.ppl0 = 6.21; % overwritten by fits keyword in reducedkeys
                    
                    % require on-sky lam/D = 9.0
                    % system lam/D = PIAAMAG * on-sky lam/D
                    % system lam/D = 6.15
                    % want on-sky control region = 9 lam/D
                    % system control region = PIAAMAG * 9 lam/D
                    
                    
                case 101 % PIAA Dan's
                    %S.Results_pn = '/proj/piaacmc/EFC/data/run101/';
                    S.Results_pn = '/proj/piaacmc/EFC/data/run101/';
                    %                     S.Rundir_pn  = 'rundir_fix/';         % always relative to Results_pn
                    %                     S.Reduced_pn = 'rundir_fix/reduced_dmarx_20210422/'; % always relative to Results_pn
                    S.Rundir_pn  = 'rundir_fix/';         % always relative to Results_pn
                    S.Reduced_pn = 'rundir_fix/reduced_dmarx_20210423/'; % always relative to Results_pn

                    
                    S.S383temp_pn= '/home/dmarx/HCIT/PIAA/hcim_testbed_run101/results/';
                    
                    S.XlimDefault = 12 * [-1 1];
                    S.YlimDefault = 12 * [-1 1];

                    S.PIAAMAG = 1.12; % should get this from config
                    S.DrawradiiDefault = S.PIAAMAG*[1.8 9.0];
                    
                    S.RminSc    = S.PIAAMAG * 1.8; % lam/D
                    S.RmaxSc    = S.PIAAMAG * 9.0;

                    % overwritten if camera image is found
                    S.NKTupper = [635]*S.NM; %[533.5, 555.5, 577.5]*S.NM;
                    S.NKTlower = [635]*S.NM; %[522.5, 544.5, 566.5]*S.NM;
                    S.NKTcenter = mean([S.NKTupper; S.NKTlower]);

                    S.ppl0 = 6.15; %4.96;

                case 102 % PIAA model 2
                    S.Results_pn = '/proj/piaacmc/EFC/data/run102/';
                    S.S383temp_pn= '/home/dmarx/HCIT/PIAA/hcim_testbed_run102/results/';
                    
                    S.XlimDefault = 12 * [-1 1];
                    S.YlimDefault = 12 * [-1 1];

                    S.PIAAMAG = 1.12; % should get this from config
                    S.DrawradiiDefault = S.PIAAMAG*[1.8 9.0];
                    
                    S.RminSc    = S.PIAAMAG * 1.8; % lam/D
                    S.RmaxSc    = S.PIAAMAG * 9.0;

                    % overwritten if camera image is found
                    S.NKTupper = [628.65]*S.NM; %[533.5, 555.5, 577.5]*S.NM;
                    S.NKTlower = [641.35]*S.NM; %[522.5, 544.5, 566.5]*S.NM;
                    S.NKTcenter = mean([S.NKTupper; S.NKTlower]);

                    
                    
                    
                case 603, % SPC_disc
                    S.Results_pn = '/home/dmarx/HCIT/SPC_disc/hcim_testbed_20170705/results/run603/';
                    S.XlimDefault = 22 * [-1 1];
                    S.YlimDefault = 22 * [-1 1];
                    
                    throughput_fn = '/home/dmarx/HCIT/SPC_disc/hcim_testbed_20170705/results/Throughput_20171003T122253_20190829.mat';
                    S.Sthpt = load(PathTranslator(throughput_fn));
                    S.Sthpt.ThptCal_fn = throughput_fn;

                    S.ppl0 = 4.01; % config_SPCdisc_20180321.py
                    S.DrawradiiDefault = [6.5 19.0];
                    
                case 604, % SPC_disc
                    S.Results_pn = '/home/dmarx/HCIT/SPC_disc/hcim_testbed_20170705/results/run604/';
                    S.XlimDefault = 22 * [-1 1];
                    S.YlimDefault = 22 * [-1 1];
                
                case 606, % MCB-SPC
                    S.Results_pn = '/home/dmarx/ln_mcb_data/EFC/SPC/run606/';
                    S.ppl0 = 6.09; % MCB SPC from config_MCB_SPC_20181015.py
                    S.XlimDefault = 12 * [-1 1];
                    S.YlimDefault = 12 * [-1 1];
                    S.DrawradiiDefault = [2.6 9.0];
                    S.DrawthetaDefault = 65*[-0.5 0.5]*CConstants.P;

                    S.RminSc    = 2.6; % lam/D
                    S.RmaxSc    = 9.0;
                    S.ThminSc   = (90 - 32.5)*CConstants.P;
                    S.ThmaxSc   = (90 + 32.5)*CConstants.P;

                    pntmp = '/home/dmarx/HCIT/MCB_SPC/hcim_testbed_run606/results/FOVThroughputMap/';
                    %ThptCal_fn = [pntmp 'fov_20181102T110148_Mjk_Thpt.mat'];
                    ThptCal_fn = [pntmp 'fov_20181102T110148_FluxHM_Thpt.mat'];
                    %S.Sthpt.fovx(:), S.Sthpt.fovy(:), S.Sthpt.thpt(:)
                    S.Sthpt = load(PathTranslator(ThptCal_fn));
                    S.Sthpt.ThptCal_fn = ThptCal_fn;

                case 607, % MCB-SPC model 2
                    S.Results_pn = '/home/dmarx/HCIT/MCB_SPC/hcim_model2_20181021/results/run607/';
                    S.ppl0 = 6.09;
                    S.XlimDefault = 12 * [-1 1];
                    S.YlimDefault = 12 * [-1 1];
                    S.DrawradiiDefault = [2.6 9.0];
                    S.DrawthetaDefault = 65*[-0.5 0.5]*CConstants.P;
                    
                    S.RminSc    = 2.6; % lam/D
                    S.RmaxSc    = 9.0;
                    S.ThminSc   = (90 - 32.5)*CConstants.P;
                    S.ThmaxSc   = (90 + 32.5)*CConstants.P;

                case 608, % MCB-SPC-IFS
                    S.Results_pn = '/home/dmarx/ln_mcb_data/EFC/SPC/run608/';

                    S.ppl0 = 6.13; % MCB SPC from config_MCB_SPC_20181015.py
                    %S.ppl0 = 3.628; % IFS value
                    % starting itnum 854, defined in ReducedKeys

                    S.XlimDefault = 12 * [-1 1];
                    S.YlimDefault = 12 * [-1 1];
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
            if exist('sOptin','var') && isa(sOptin,'struct'),
                ffields = fieldnames(sOptin);
                for ii = 1:length(ffields),
                    S.(ffields{ii}) = sOptin.(ffields{ii});
                end
            end
     
            % build paths and filenames for the data
            % basename, such as 'run603it00000.fits'
            s_bn = ['run' num2str(S.runnum,'%03d') 'it' num2str(iter,'%05d') '.fits'];
            S.Reduced_fn = PathTranslator(cell2mat(join({S.Results_pn, S.Reduced_pn, s_bn}, '/')));
            S.Rundir_fn = PathTranslator(cell2mat(join({S.Results_pn, S.Rundir_pn, s_bn}, '/')));

            % 'local' paths are either the same (linux cluster)
            %         or local PC drive. If PC, we will copy data
            %         locally to save unzip and read time
            if ispc,
                if ~exist([S.PCtemp_pn S.Reduced_pn],'dir')
                    mkdir([S.PCtemp_pn S.Reduced_pn])
                end
                sReduced_local_fn = fullfile(S.PCtemp_pn, S.Reduced_pn, s_bn);
                sRundir_local_fn  = fullfile(S.PCtemp_pn, S.Rundir_pn, s_bn);
            else
                if ~exist(fullfile(S.S383temp_pn, S.Reduced_pn),'dir')
                    mkdir(fullfile(S.S383temp_pn, S.Reduced_pn))
                end
                sReduced_local_fn = fullfile(S.S383temp_pn, S.Reduced_pn, s_bn);
                sRundir_local_fn  = fullfile(S.S383temp_pn, S.Rundir_pn, s_bn);
            end
            
            % unzip the fits files if necessary
            % first check if 'local' fits files exists,
            try,
                if ~exist(sReduced_local_fn, 'file'),
                    sReduced_gz = [sReduced_local_fn '.gz'];
                    if ~exist(sReduced_gz, 'file') %&& ispc,
                        copyfile([S.Reduced_fn '.gz'], sReduced_gz);
                    end
                    gunzip(sReduced_gz);
                    % now reduced data should be unzipped, and
                    % S.Reduced_fn = the 'local' unzipped data file
                end
                S.Reduced_fn = sReduced_local_fn;
            catch ME
                disp('File Error:');
                disp(ME.identifier);
                disp([S.Reduced_fn '.gz']);
                % return empty instance
                return
            end 
            % repeat local logic for rundir
            if ~exist(sRundir_local_fn, 'file'),
                sRundir_gz = [sRundir_local_fn '.gz'];
                if ~exist(sRundir_gz, 'file') %&& ispc,
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
            
            % useful parameters from ReducedKeys
            ppl0tmp = FitsGetKeywordVal(S.ReducedKeys,'ppl0');
            if ~isempty(ppl0tmp), S.ppl0 = ppl0tmp; end
            
            % bscan list
            if ~isempty(S.ReducedKeys)
                ireg = 0;
                while true
                    btmp = FitsGetKeywordVal(S.ReducedKeys, ['BSCAN' num2str(ireg,'%03d')]);
                    if isempty(btmp), break, end
                    S.bscan(ireg+1) = btmp;
                    ireg = ireg + 1;
                end
                S.Nbscan = length(S.bscan);
            end
            S.betamin = FitsGetKeywordVal(S.ReducedKeys, 'BMIN');
            
            
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
            S.imgindex = 1:S.NumImProbe:S.NofW*S.NumImProbe; % slice of ImCube that is unprobed for each wvl
            
            % get actual wavelengths from the camera fits files
            for iwv = 1:S.NofW,
               fn = FitsGetKeywordVal(S.ImKeys, ['C' num2str(iwv-1) 'P0J0']) ;
               if isempty(fn),
                   warning(['no header key ' ['C' num2str(iwv-1) 'P0J0'] ]);
                   break;
               end
               
               [pntmp, fntmp, ext] = fileparts(fn);
               fn = [pntmp '/' fntmp '.fits'];
               if isempty(fn) || ~exist(PathTranslator(fn),'file'),
                   warning(['cannot open ' fn]);
                   %continue
                   % hack for now
                   %fn = ['/home/dmarx/ln_mcb_data/IFS/images/' fn];
               else
                   finfo = fitsinfo(PathTranslator(fn));
                   if ~isempty(FitsGetKeywordVal(finfo.PrimaryData.Keywords,'NKTLOWER')),
                       S.NKTlower(iwv) = FitsGetKeywordVal(finfo.PrimaryData.Keywords,'NKTLOWER')*S.NM;
                       S.NKTupper(iwv) = FitsGetKeywordVal(finfo.PrimaryData.Keywords,'NKTUPPER')*S.NM;
                       S.NKTcenter(iwv) = mean([S.NKTlower(iwv) S.NKTupper(iwv)]);
                   end
               end
            end
            
            S.lambda = S.NKTcenter;

            % get timestamp from raw image file
            fntmp = FitsGetKeywordVal(S.ImKeys, ['C' num2str(iwv-1) 'P0J0']);
            if ~isempty(fntmp)
                dtmp = dir(PathTranslator(fn));
                S.timestamp = datetime(dtmp.date);
            end
    
            % if requested in varargin, run methods right away
            for icom = 1:length(varargin),
                if any(strcmpi(varargin{icom}, methods(S)))
                    S.(varargin{icom})(varargin{icom+1:end});
                end
            end
                
        end % function CRunData
        
        function Contrast = GetContrast(S, varargin)
            % % options
            % rminsc = CheckOption('RminSc', S.RminSc, varargin{:});
            % rmaxsc = CheckOption('RmaxSc', S.RmaxSc, varargin{:});
            % yminsc = CheckOption('YminSc', -inf, varargin{:});
            % ymaxsc = CheckOption('YmaxSc', inf, varargin{:});
            % xminsc = CheckOption('XminSc', -inf, varargin{:});
            % xmaxsc = CheckOption('XmaxSc', inf, varargin{:});
            % bMaskScUse= CheckOption('bMaskSc', S.bMaskSc, varargin{:});
            % bDisplay = CheckOption('display', true, varargin{:});
            %
            % rminsc, rmaxsc, etc are applied to bMaskSc
            
            % initialize empty return struct
            Contrast = struct( ...
                'cntl_lam', [] ...   % unprobed image control region NI
                ,'score_lam', [] ... % unprobed image score region NI
                ,'contr_lam', [] ... % contrast score region
                ,'mean', [] ...      % mean(contr_lam)
                ,'inco_lam', [] ... % score region
                ,'co_lam', [] ...   % score region
                ,'inco_lam_NI', [] ... % score region
                ,'co_lam_NI', [] ...   % score region
                ,'inco_mean', [] ...   % score region
                ,'co_mean', [] ...     % score region
                );
            
            if isempty(S.bMask),
                S.ReadMaskCube;
            end
            
            if isempty(S.ImCube),
                S.ReadImageCube;
            end
            
            if isempty(S.CohInt),
                S.ReadReducedCube;
            end

            if isempty(S.NofW),
                % this is an empty instance
                return
            end
            
            % options
            rminsc = CheckOption('RminSc', S.RminSc, varargin{:});
            rmaxsc = CheckOption('RmaxSc', S.RmaxSc, varargin{:});
            yminsc = CheckOption('YminSc', S.YminSc, varargin{:});
            ymaxsc = CheckOption('YmaxSc', S.YmaxSc, varargin{:});
            xminsc = CheckOption('XminSc', S.XminSc, varargin{:});
            xmaxsc = CheckOption('XmaxSc', S.XmaxSc, varargin{:});
            bMaskScUse= CheckOption('bMaskSc', S.bMaskSc, varargin{:});
            bDisplay = CheckOption('display', true, varargin{:});
            
            % coordinate system in back-end lam/D
            [x, y, X, Y, R, T] = CreateGrid(bMaskScUse, 1./S.ppl0);

            % FOV score mask for calculating contrast
            % options can override object settings
            bMaskScUse = bMaskScUse & (R >= rminsc) & (R <= rmaxsc) & (Y >= yminsc) & (Y <= ymaxsc) & (X >= xminsc) & (X <= xmaxsc);

            % resample throughput data to pixels in scoring region
            if ~isempty(S.Sthpt),

                %                 if S.runnum == 603, % SPC disc
                %                     error('SPC disc throughput data format needs to be reworked');
                %                     %                     dxyfudge = 0.9;
                %                     %                     rthpt = mean([Sth.src_dy(:).*dxyfudge./Sth.mmperlamD Sth.src_dx(:).*dxyfudge./Sth.mmperlamD],2);
                %                     %                     thpt = mean([Sth.thpt_dy(:) Sth.thpt_dx(:)],2);
                %                     %                     % normalize to clear region: lam/D > 8.5 & lam/D < 16.2
                %                     %                     thpt_norm = mean(thpt(abs(rthpt) >= 8.5 & abs(rthpt) < 16.2));
                %                     %                     thpt = thpt./ thpt_norm;
                %                     %                     % rthpt, thpt is now a lookup table to throughput for
                %                     %                     % cross-section
                %                     %                     % make a matrix
                %                     %                     Thpt = zeros(size(bMaskSc));
                %                     %                     Thpt(bMaskSc) = interp1(rthpt, thpt, R(bMaskSc));
                %                 end
                if isfield(S.Sthpt,'fovx') && isfield(S.Sthpt,'fovy'),
                    % most cases
                    Finterp = scatteredInterpolant(S.Sthpt.fovx(:), S.Sthpt.fovy(:), S.Sthpt.thpt(:));
                    S.Sthpt.ThptSc = ones(size(bMaskScUse)); % avoid divide by zero
                    S.Sthpt.ThptSc(bMaskScUse) = Finterp(X(bMaskScUse), Y(bMaskScUse));
                elseif isfield(S.Sthpt,'fovr'),
                    % SPC Disc run 603
                    S.Sthpt.ThptSc = ones(size(bMaskScUse)); % avoid divide by zero
                    S.Sthpt.ThptSc(bMaskScUse) = interp1(S.Sthpt.fovr, S.Sthpt.thpt, R(bMaskScUse));                    
                end

            else % S.Sthpt is empty
                %warning('no throughput data, contrast is normalized intensity');
                S.Sthpt.ThptSc = ones(size(bMaskScUse));
            end
            
            % note: don't like nonzeros(...) because pixels within the dark
            % hole that happen to have zero measured intesnity do not
            % contribute to the mean.
            Contrast.cntl_lam = zeros(1,S.NofW);
            for iwv = 1:S.NofW,
                % NI total control region
                Contrast.cntl_lam(iwv) = mean(nonzeros( S.ImCube(:,:,S.imgindex(iwv)).*S.bMask ));
                % NI total score region
                Contrast.score_lam(iwv) = mean(S.ImCubeUnProb{iwv}(bMaskScUse)); 
                % Contrast score region, Modulated & Unmodulated as
                % contrast
                S.ImCubeContrast{iwv} = S.ImCubeUnProb{iwv}./S.Sthpt.ThptSc;
                S.ImCubeCohContrast{iwv} = S.CohInt{iwv}./S.Sthpt.ThptSc;
                S.ImCubeIncContrast{iwv} = S.IncInt{iwv}./S.Sthpt.ThptSc;
                Contrast.contr_lam(iwv) = mean(S.ImCubeContrast{iwv}(bMaskScUse)./S.Sthpt.ThptSc(bMaskScUse));
                
            end
            Contrast.mean = mean(Contrast.contr_lam);

            % mean coherent and incoherent contrast
            
            % scoring (or user-defined) region
            Contrast.inco_lam = zeros(1,S.NofW);
            Contrast.co_lam = zeros(1,S.NofW);
            for iwv = 1:S.Nlamcorr,
                
                % total incoherent, use IncIntEst:
                % pixels where inc int < 0, fixed to = eps
                Contrast.inco_lam_NI(iwv) = mean(S.IncIntEst{iwv}(bMaskScUse));
                Contrast.inco_lam(iwv) = mean(S.IncIntEst{iwv}(bMaskScUse)./S.Sthpt.ThptSc(bMaskScUse));
                
                % Mix incoherent (those pixels where sigtbw < 0, Eest = 0)
                
                % coherent, use S.CohIntEst, pixels where inc int < 0, coh int is fixed to = unprobed
                Contrast.co_lam_NI(iwv) = mean(S.CohIntEst{iwv}(bMaskScUse));
                Contrast.co_lam(iwv)   = mean(S.CohIntEst{iwv}(bMaskScUse)./S.Sthpt.ThptSc(bMaskScUse));
            end
            Contrast.inco_mean = mean(Contrast.inco_lam);
            Contrast.co_mean = mean(Contrast.co_lam);

            % radial plot of contrast
            if bDisplay,
                [hfigrad, harad] = S.DisplayRadialPlot(S.ImCubeContrast, ...
                    'ylabel', 'Contrast', 'dispradlim', [0 max(R(bMaskScUse))]);
            end
            
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
                %fprintf('no reduced data\n');
                return
            end
            
            if isempty(S.ImCube),
                S.ReadImageCube;
            end

            
            try,
                RedData = fitsread(S.Reduced_fn); % primary hdu
            catch ME
                disp('File Error:');
                disp(ME.identifier);
                disp([S.Reduced_fn]);
                return                
            end

            
            [nr, nc, nslices] = size(RedData);
            % each images slice is nr x nc
            
            % initialize some data cubes
            S.bPampzero    = false(S.Nlamcorr, nr, nc);
            S.bMask_badpix = false(S.Nlamcorr, nr, nc);
            
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

            if isempty(S.bMask),  S.ReadMaskCube; end
            for iwl = 1:S.Nlamcorr,

                % incoherent image
                S.IncInt{iwl} = RedData(:,:,0*S.Nlamcorr+iwl);

                % probe fields
                for ipr = 1:S.Nppair,
                    S.ProbeAmp{iwl,ipr} = RedData(:,:,ipr*S.Nlamcorr+iwl);
                end

                % coherent intensity
                S.CohInt{iwl} = RedData(:,:,(1+S.Nppair)*S.Nlamcorr+iwl).^2 + RedData(:,:,(1+S.Nppair+1)*S.Nlamcorr+iwl).^2;

                % coherent complex estimated field and model field
                S.E_t(iwl,:,:)	= (RedData(:,:,(1+S.Nppair)*S.Nlamcorr+iwl)+1i*RedData(:,:,(1+S.Nppair+1)*S.Nlamcorr+iwl));
                S.E_m(iwl,:,:) 	= (RedData(:,:,(3+S.Nppair)*S.Nlamcorr+iwl)+1i*RedData(:,:,(3+S.Nppair+1)*S.Nlamcorr+iwl)); % starting efield

                % dedall optimal beta,
                S.dE_optimal(iwl,:,:) = (RedData(:,:,(5+S.Nppair)*S.Nlamcorr+iwl)+1i*RedData(:,:,(5+S.Nppair+1)*S.Nlamcorr+iwl));
                % dedall each bscan
                for ireg = 1:S.Nbscan,
                    nnn = 2*(ireg-1);
                    S.dE_bscan(iwl,:,:) = (RedData(:,:,(7+nnn+S.Nppair)*S.Nlamcorr+iwl)+1i*RedData(:,:,(7+nnn+S.Nppair+1)*S.Nlamcorr+iwl));
                end

                % pampzero is logical, imag and real are the same logical
                % pampzero = true for pixels where: eest == 0 OR iinc < -()
                %    OR eestcond < eestcondlim
                nnn = 2*S.Nbscan;
                S.bPampzero(iwl,:,:) = logical(RedData(:,:,(7+nnn+S.Nppair)*S.Nlamcorr+iwl));

                % recreate bPampzero, note: this might differ from the
                % pampzero used in EFC
                % don't use pixels where, abs(CohInt)==0, IncInt < 0
                % pampzero = false for good pixels, mdMask is control
                % region, bMask is logical(mdMask)
                S.bMask_badpix(iwl,:,:) = ~S.bMask | abs(S.CohInt{iwl}) == 0 | S.IncInt{iwl} < 0;
                
                %%%% NOTE: NormIntensity_ does not account for imwt, use
                %%%% GetContrast() for more options
                % mean coherent and incoherent contrast
                S.NormIntensity_inco(iwl) = mean(S.IncInt{iwl}(~S.bMask_badpix(iwl,:,:)));
                S.NormIntensity_co(iwl)   = mean(S.CohInt{iwl}(~S.bMask_badpix(iwl,:,:)));

                % if require all probes is true, pixels where abs(CohInt) == 0,
                % are pixels where probe amp < 0 for at least one probe,
                % image at this pixel is considered all incoherent.
                % Separate into two components
                bProbNeg = S.bMask .* ~( abs(S.CohInt{iwl}) > 0 );
                S.IncIntMix{iwl} = (~bProbNeg).*S.IncInt{iwl}; % inc int = 0 at these points
                
                % if probe estimates coherent intensity is larger than
                % unprobed intensity, inc int is given intensity < 0
                % have option to just set inc int to zero at these poitns
                S.IncIntEst{iwl} = S.IncInt{iwl};
                S.IncIntEst{iwl}(S.IncInt{iwl} <= 0) = eps;
                % and then where inc int is zero, coh int = unprobed
                S.CohIntEst{iwl} = S.CohInt{iwl};
                S.CohIntEst{iwl}(S.IncInt{iwl} <= 0) = S.ImCubeUnProb{iwl}(S.IncInt{iwl} <= 0);
                                       
            end % for each wl
            
            % whole band mean
            S.CohIntFullBand = mean(cat(3, S.CohIntEst{:}), 3);
            S.IncIntEstFullBand = mean(cat(3, S.IncIntEst{:}), 3);
            S.IncIntFullBand = mean(cat(3, S.IncInt{:}), 3);
            
            
        end % ReadReducedCube
        
        function S = ReadSVD(S)
            % if svd spectrum was saved
            fntmp = PathTranslator([S.Results_pn, S.Reduced_pn, ['svd_it' num2str(S.iter, '%05d') '.mat']]);
            
            try
                if exist(fntmp, 'file')
                    S.svd_fn = fntmp;
                    S.sSVD = load(S.svd_fn);
                end
                
                S.sSVD.s2norm = double((1./S.sSVD.mjsvs).*(S.sSVD.s.^2));
                S.sSVD.SpecIntensity = double(abs(S.sSVD.U'*S.sSVD.rhs).^2);

            catch ME
                disp('File Error:');
                disp(ME.identifier);
                disp([S.Reduced_fn]);
                
            end
                        
        end % ReadSVD
        
        function S = ReadMaskCube(S, varargin)
            
            % 1st slice is mask of control region (md)
            % next nlamcorr * 2 slices are real and imag occulter transmission profile

            try,
                MaskCubeData = fitsread(S.Reduced_fn,'image',1);
            catch ME
                disp('File Error:');
                disp(ME.identifier);
                disp([S.Reduced_fn]);
                return
            end
            
            % 1st slice is mask of control region, i.e. model.band[ilam0].md
            S.mdMask = squeeze(MaskCubeData(:,:,1));
            S.bMask  = S.mdMask > 0;
            
            % let's create score mask as well
            % default values from the class
            %val = CheckOption(sOpt, valDefault, varargin)
            RminSc = CheckOption('RminSc', S.RminSc, varargin{:});
            RmaxSc = CheckOption('RmaxSc', S.RmaxSc, varargin{:});
            TminSc = CheckOption('TminSc', S.ThminSc, varargin{:});    
            TmaxSc = CheckOption('TmaxSc', S.ThmaxSc, varargin{:});                
            YminSc = CheckOption('YminSc', S.YminSc, varargin{:});
            YmaxSc = CheckOption('YmaxSc', S.YmaxSc, varargin{:});
            XminSc = CheckOption('XminSc', S.XminSc, varargin{:});
            XmaxSc = CheckOption('XmaxSc', S.XmaxSc, varargin{:});
            
            [x, y, X, Y, R, T] = CreateGrid(S.bMask, 1./S.ppl0);
            S.bMaskSc = S.bMask & (R >= RminSc & R <= RmaxSc & Y >= YminSc & Y <= YmaxSc & X >= XminSc & X <= XmaxSc);
            if ~any(S.bMaskSc(:)),
                warning('score region mask is empty!');
            end            
            
            if ~isempty(TminSc) && ~isempty(TmaxSc),
                S.bMaskSc = S.bMaskSc & ...
                    ( (T >= TminSc & T <= TmaxSc) | (T >= mod2pi(TminSc + pi) & T <= mod2pi(TmaxSc + pi)) );
                
                %figure, imageschcit(x, y, bMaskSc), title('Check Scoring Region Mask')
            end
            
            
            
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
            
            % finfo = fitsinfo(S.Reduced_fn);
            % ProbeKwds = finfo.Image(2).Keywords; % nothing interesting
            ProbeData = fitsread(S.Reduced_fn,'image',2);

            for iwl = 1:S.Nlamcorr,
                for ip = 1:S.Nppair,
                    
                    S.ProbeModel{iwl, ip} = squeeze(...
                        ProbeData(:,:,0*S.Nlamcorr*S.Nppair+(ip-1)*S.Nlamcorr+iwl) ...
                        .* exp(1i*...
                        ProbeData(:,:,1*S.Nlamcorr*S.Nppair+(ip-1)*S.Nlamcorr+iwl)) ...
                        );
                    
                    S.ProbeMeasAmp{iwl, ip} = real(sqrt(squeeze( ProbeData(:,:,2*S.Nlamcorr*S.Nppair + (ip-1)*S.Nlamcorr+ iwl) )));
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
            
            if isempty(S.bMask), S.ReadMaskCube; end
            
            try,
                S.ImCube = fitsread(S.Rundir_fn);
            catch ME
                disp('File Error:');
                disp(ME.identifier);
                disp([S.Rundir_fn]);
                return                
            end
            
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
                
                S.NormIntensity_total(iwv) = mean(S.ImCubeUnProb{iwv}(S.bMask)); % bMask = logical(mdMask)
                
            end
            
            % Useful Derived Parameters                       
            % mean un probed over the band
            S.ImCubeUnProbFullBand = mean(cat(3, S.ImCubeUnProb{:}), 3);
            
        end % ReadImageCube
        
        function S = ReadDMvCube(S)
            
            try,
                finfo = fitsinfo(S.Rundir_fn);
            catch ME
                disp('File Error:');
                disp(ME.identifier);
                disp(S.Rundir_fn);
                return                
            end

            % Primary hdu is unprobed images
            % one Image hdu per DM
            S.Ndm = length(finfo.Image);
            
            for dmnum = 1:S.Ndm,
                S.DMvCube{dmnum} = fitsread(S.Rundir_fn,'image',dmnum);
            end
            
        end % ReadDMv

        function [hfig, hax, sImageCubeData] = DisplayImCubeImage(S, varargin)
            % [hfig, him, Im] = DisplayImCubeImage(S)
            % [hfig, him, Im] = DisplayImCubeImage(S, imnum)
            % 
            % raw camera images, adjusted for photometry and dark
            % these are the images from rundir/*.fits.gz
            %
            % display all, imnum, in an ImageCube
            % 
            % Options:
            %   CheckOption('scale', 'log', varargin{:}); % 'linear', 'log'
            %   CheckOption('imnum', [], varargin{:});
            %   CheckOption('fTitleStr', @(isl) ['Iter #' num2str(S.iter) ', Image #' num2str(imnumlist(isl))]);
            %
            % and varargin passed to ImageCube
            
            logscale = CheckOption('scale', 'log', varargin{:}); % 'linear', 'log'
            imnumlist = CheckOption('imnum', [], varargin{:});
            
            if isempty(S.ImCube),
                S.ReadImageCube;
            end
            
            % make image cube
            imcube = shiftdim(S.ImCube,2);
            [Nsl, nr, nc] = size(imcube);
            if isempty(imnumlist),
                imnumlist = 1:Nsl;
            end
            imcube = imcube(imnumlist,:,:);
            fTitleStr = CheckOption('fTitleStr', @(isl) ['Iter #' num2str(S.iter) ', Image #' num2str(imnumlist(isl))], varargin{:});
        
            % scale
            switch logscale
                case 'linear'
                    % do nothing
                case 'log'
                    imcube = log10((imcube>0).*imcube);
                    clim = max(imcube(:)) + [-3 0];
                otherwise
                    error(['unknown scale ' logscale]);
            end
            
            figure;
            [hfig, hax, sImageCubeData] = ImageCube(imcube, imnumlist, 'fTitleStr', fTitleStr, varargin{:});
            set(hax,'clim',clim)
            colorbar
            
        end % DisplayImCubeImage
        
        function [hfig, hax, x, y] = DisplaySingleImage(S, Im, varargin)
            % S.DisplaySingleImage(Im, options)
            %
            % generic display an image using the standard format for this
            % run, etc.
            %
            % default options and set requested options
            %             bPlotLog = CheckOption('bLog', false, varargin{:});
            %             sTitle = CheckOption('title', '', varargin{:});
            %             drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            %             drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            %             climopt = CheckOption('clim', [], varargin{:});
            %             xcoff = CheckOption('xcoff', 0, varargin{:}); % star offset, pixels
            %             ycoff = CheckOption('ycoff', 0, varargin{:}); % star offset, pixels
            
            % default options and set requested options
            %  val = CheckOption(sOpt, valDefault, varargin)
            bPlotLog = CheckOption('bLog', false, varargin{:});
            sTitle = CheckOption('title', '', varargin{:});
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            climopt = CheckOption('clim', [], varargin{:});
            xcoff = CheckOption('xcoff', 0, varargin{:}); % star offset, pixels
            ycoff = CheckOption('ycoff', 0, varargin{:}); % star offset, pixels
            hax = CheckOption('hax', [], varargin{:});
            hfig = CheckOption('hfig', [], varargin{:});
            
            [x, y] = CreateGrid(Im, 1./S.ppl0);
            x = x - xcoff/S.ppl0;
            y = y - ycoff/S.ppl0;
            
            if isempty(hfig) && isempty(hax),
                hfig = figure;
            end
            if ~isempty(hfig)
                figure(hfig)
            end
            if ~isempty(hax),
                axes(hax)
            end
            
            if bPlotLog,
                imageschcit(x, y, log10(abs(Im))), axis image,
                colorbartitle('log_{10} Norm Intensity')
            else
                imageschcit(x, y, Im), axis image,
                colorbartitle('Norm Intensity (linear)')
            end
            title(sTitle)
            
            hax = gca;
            
            set(hax,'xlim',xlim,'ylim',ylim)
            xlabel('\lambda/D'), ylabel('\lambda/D')
            
            if ~isempty(climopt)
                set(hax,'clim',climopt)
            end
            
            % overlay circles if requested
            DrawCircles(hax, drawRadii);
            DrawThetaLines(hax, drawTheta, drawRadii);

        end % DisplaySingleImage

        function [sHandles, rplot, IntRad] = DisplayRadialPlot(S, ImCube, varargin)
            % [sHandles, rplot, IntRad] = DisplayRadialPlot(S, ImCube, varargin)
            % generic routine for radial plot of intensity or contrast
            % ImCube is cell array (1 x NofW)
            % e.g. ImCubeUnProbe, ImCubeUnProbe/Thrpt
            %      ImCubeCohInt, ImCubeIncInt
            %
            % output rplot, IntRad = radius, radial intensity data
            %            
            %             Nr = CheckOption('nr', ceil(min([128 length(R)/4])), varargin{:}); % # of radial sample pts
            %             dispRadlim = CheckOption('dispradlim', [0 max(S.DrawradiiDefault)], varargin{:});
            %             drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            %             bMaskUse = CheckOption('bMask', S.bMask, varargin{:});
            %             strYlabel = CheckOption('ylabel', 'Average Normalized Intensity', varargin{:});
            %             plotRequired = CheckOption('plotrequired', [], varargin{:}); % [r(:) contrast(:)]
            %             iplot = CheckOption('iplot', 1:length(ImCube), varargin{:});
            %             legstr = CheckOption('legstr', [], varargin{:});
            %             strTitle = CheckOption('title', ['Iter #' num2str(S.iter)], varargin{:});
            %             bPlotMean = CheckOption('plotmean', true, varargin{:});
            %             haxuse = CheckOption('hax', [], varargin{:});

            [x, y, X, Y, R] = CreateGrid(ImCube{1}, 1./(S.ppl0*S.PIAAMAG));

            Nr = CheckOption('nr', ceil(min([128 length(R)/4])), varargin{:}); % # of radial sample pts
            dispRadlim = CheckOption('dispradlim', [0 max(S.DrawradiiDefault)], varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            bMaskUse = CheckOption('bMask', S.bMask, varargin{:});
            strYlabel = CheckOption('ylabel', 'Average Normalized Intensity', varargin{:});
            plotRequired = CheckOption('plotrequired', [], varargin{:}); % [r(:) contrast(:)]
            iplot = CheckOption('iplot', 1:length(ImCube), varargin{:});
            legstr = CheckOption('legstr', [], varargin{:});
            strTitle = CheckOption('title', ['Iter #' num2str(S.iter)], varargin{:});
            bPlotMean = CheckOption('plotmean', true, varargin{:});
            haxuse = CheckOption('hax', [], varargin{:});
            
            re = linspace(dispRadlim(1), dispRadlim(2), Nr+1)';
            IntRad = cell(length(iplot),1);
            legstrwv = cell(1,length(iplot));
            for iiwv = 1:length(iplot),
                iwv = iplot(iiwv);
                Itmp = zeros(Nr,1);
                for ir = 1:Nr,
                    % including S.bMask applies theta (bowtie) limits
                    Itmp(ir) = mean(ImCube{iwv}(R > re(ir) & R <= re(ir+1) & bMaskUse));
                end % for ir
                IntRad{iiwv} = Itmp;
                %legstrwv{iiwv} = [num2str(S.NKTcenter(iwv)/S.NM,'%.1f') 'nm'];
            end % for iwv
            rplot = mean([re(1:end-1) re(2:end)],2); % radii midway between edges
            
            if isempty(legstr), legstr = legstrwv; end            
            
            if ~isempty(haxuse),
                axes(haxuse);
                hfig = gcf;
            else,
                hfig = figure;
            end 
            hl = semilogy(rplot, ([IntRad{:}]>0).*[IntRad{:}]);
            ha = gca;
            hold on
            
            % add plot of mean
            if bPlotMean,
                hl(end+1) = semilogy(rplot, mean([IntRad{:}],2), '-k');
                legstr{end+1} = 'Mean';
                set(hl(end),'LineWidth',2);
            end
            
            % add plot of contrast requirement, if provided
            if ~isempty(plotRequired),
                hl(end+1) = semilogy(plotRequired(:,1), plotRequired(:,2), '--r');
                set(hl(end),'LineWidth',2);
                legstr{end+1} = 'Requirement';
            end %
            
            grid on
            
            if ~isempty(drawRadii),
                ylim = get(gca,'ylim');
                hold on
                for irad = 1:length(drawRadii),
                    plot(drawRadii(irad)*[1 1], ylim, '--r')
                    %legstr{end+1} = [num2str(drawRadii(irad),'%.1f')];
                end
                hold off
            end

            xlabel('Radius (\lambda/D)')
            ylabel(strYlabel)
            hleg = legend(legstr{:}); %, 'location','north');
            title(strTitle)

            % return all the handles
            sHandles = struct(...
                'hfig', hfig ...
                ,'hax', ha ...
                ,'hline', hl ...
                ,'hleg', hleg ...
                );
            
        end % DisplayRadialPlot
        
        function [hfig, ha] = DisplayImCubeUnProb(S, varargin)
            % [hfig, ha] = DisplayImCubeUnProb(S, ...)
            % default options and set requested options
            %             bPlotLog = CheckOption('bLog', false, varargin{:});
            %             drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            %             drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            %             drawylimlines = CheckOption('drawylimlines', [], varargin{:})
            %             climopt = CheckOption('clim', [], varargin{:});
            %             ilam = CheckOption('ilam', 1:S.NofW, varargin{:});
            %             haxuse = CheckOption('hax', [], varargin{:});
            %             xlim = CheckOption('xlim', [], varargin{:});
            %             ylim = CheckOption('ylim', [], varargin{:});
            %
            
            if isempty(S.ImCubeUnProb),
                try
                    S.ReadImageCube;
                catch ME
                    disp(ME.message);
                    hfig = [];
                    ha = [];
                    return
                end
            end
            %             if isempty(S.bMask),
            %                 S.ReadMaskCube;
            %             end
            
            % default options and set requested options
            %  val = CheckOption(sOpt, valDefault, varargin)
            bPlotLog = CheckOption('bLog', false, varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            drawYlimLines = CheckOption('drawylimlines', [], varargin{:});
            climopt = CheckOption('clim', [], varargin{:});
            ilam = CheckOption('ilam', 1:S.NofW, varargin{:});
            haxuse = CheckOption('hax', [], varargin{:});
            bPlotRadialIntensity = CheckOption('DisplayRadialIntensity', false, varargin{:});
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            
            [x, y, X, Y, R] = CreateGrid(S.ImCubeUnProb{1}, 1./(S.ppl0*S.PIAAMAG));
            
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
                    %imageschcit(x/S.PIAAMAG, y/S.PIAAMAG, log10(abs(S.ImCubeUnProb{iwvpl}))), axis image,
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
            DrawYlimLines(ha, drawYlimLines, drawRadii);
            % DrawYlimLines                        CheckOption('drawylimlines', [], varargin{:})
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
            if bPlotRadialIntensity,
                [hfigrad, harad] = S.DisplayRadialPlot(S.ImCubeUnProb, varargin{:});
            end
            
        end % DisplayImCubeUnProb(S)
        
        function [hfig, ha] = DisplayImCubeContrast(S, varargin)
            % [hfig, ha] = DisplayImCubeContrast(S, ...)
            % default options and set requested options
            %             bPlotLog = CheckOption('bLog', false, varargin{:});
            %             xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            %             ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            %             drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            %             drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            %             climopt = CheckOption('clim', [], varargin{:});
            %             ilam = CheckOption('ilam', 1:S.NofW, varargin{:});
            %             haxuse = CheckOption('hax', [], varargin{:});
            %
            
            if isempty(S.ImCubeContrast),
                S.GetContrast;
            end
            if isempty(S.bMask),
                S.ReadMaskCube;
            end
            
            % default options and set requested options
            %  val = CheckOption(sOpt, valDefault, varargin)
            bPlotLog = CheckOption('bLog', true, varargin{:});
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});            
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            climopt = CheckOption('clim', [], varargin{:});
            ilam = CheckOption('ilam', 1:S.NofW, varargin{:}); % which wavelengths to plot
            haxuse = CheckOption('hax', [], varargin{:});
            
            [x, y, X, Y, R] = CreateGrid(S.ImCubeContrast{1}, 1./S.ppl0);
            
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
                    imageschcit(x, y, log10(abs(S.ImCubeContrast{iwvpl}.*S.bMaskSc))), axis image,
                    colorbartitle('log_{10} Contrast')
                else
                    imageschcit(x, y, S.ImCubeContrast{iwvpl}.*S.bMaskSc), axis image,
                    colorbartitle('Contrast')
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

            % quick hack to plot total mean in one plot
            figure,
            ImMean = zeros(size(S.ImCubeContrast{1}));
            for iwv = 1:S.NofW,
                ImMean = ImMean + S.ImCubeContrast{iwv};
            end
            ImMean = ImMean ./ S.NofW;
            imageschcit(x, y, log10(abs(ImMean.*S.bMaskSc))), axis image,
            set(gca,'xlim',xlim,'ylim',ylim)
            xlabel('\lambda/D'), ylabel('\lambda/D')
            colorbartitle('log_{10} Contrast')
            DrawCircles(gca, drawRadii);
            DrawThetaLines(gca, drawTheta, drawRadii);
            set(gca,'clim',get(ha(1),'clim'))
            title('Mean Contrast Full Band')
            
        end % DisplayImCubeContrast(S)

        function [hfig, ha] = DisplayImCubeSigProbes(S, varargin)
            % [hfig, ha] = DisplayImCubeSigProbes(S, varargin)
            %
            % S.ImCubeSigProb{iwv,ipr} = 0.5*(ImPrPlus + ImPrMinus) - S.ImCubeUnProb{iwv};
            %
            % options:
            %             bPlotLog = CheckOption('bLog', false, varargin{:});
            %             xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            %             ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            %             drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            %             drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            %             climopt = CheckOption('clim', [], varargin{:});
            %             ilam = CheckOption('ilam', 1:S.NofW, varargin{:});
            %             ipro = CheckOption('iprobe', 1:S.Nppair, varargin{:});
            %             haxuse = CheckOption('hax', [], varargin{:});
            
            if isempty(S.ImCube),
                S.ReadImageCube;
            end
            
            % default options and set requested options
            %  val = CheckOption(sOpt, valDefault, varargin)
            bPlotLog = CheckOption('bLog', false, varargin{:});
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});            
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            drawTheta = CheckOption('drawtheta', S.DrawthetaDefault, varargin{:});
            climopt = CheckOption('clim', [], varargin{:});
            ilam = CheckOption('ilam', 1:S.NofW, varargin{:});
            ipro = CheckOption('iprobe', 1:S.Nppair, varargin{:});
            haxuse = CheckOption('hax', [], varargin{:});
            
            [x, y, X, Y, R] = CreateGrid(S.ImCubeSigProb{1,1}, 1./S.ppl0);
            
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
            %             CheckOption('bLog', true, varargin{:});
            %             CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            %             CheckOption('drawylimlines', [], varargin{:});
            %             CheckOption('clim', [], varargin{:});
            %             CheckOption('hax', [], varargin{:}); % put image on this axes
            %             CheckOption('type', 'est', varargin{:});
            %             CheckOption('xlim', [], varargin{:});
            %             CheckOption('ylim', [], varargin{:});
            
            if isempty(S.IncInt),
                S.ReadReducedCube;
            end
            
            % options:
            % %  val = CheckOption(sOpt, valDefault, varargin)
            bLog = CheckOption('bLog', true, varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            drawYlimLines = CheckOption('drawylimlines', [], varargin{:});
            clim = CheckOption('clim', [], varargin{:});
            haxuse = CheckOption('hax', [], varargin{:}); % put image on this axes
            IncIntType = CheckOption('type', 'est', varargin{:});
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            
            %%%% end options

            switch lower(IncIntType),
                case 'normal',
                    plIncInt = S.IncInt;
                case 'est',
                    plIncInt = S.IncIntEst; % IncInt(IncInt < 0) = eps
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
            
            if isempty(clim),
                clim = pClim(pFun([plIncInt{:}]));
            end
            
            [x, y] = CreateGrid(plIncInt{1}, 1./(S.ppl0*S.PIAAMAG));
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
            DrawYlimLines(hax, drawYlimLines, drawRadii);
            
        end % DisplayIncInt
                
        function [hfig, hax] = DisplayProbeAmp(S, varargin)
            if isempty(S.ProbeAmp),
                S.ReadReducedCube;
            end
            
            % check options
            clim = CheckOption('clim', [], varargin{:});
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            
            % grid in lam/D
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
                    set(gca,'xlim', xlim, 'ylim', ylim)
                    title(['Iteration # ' num2str(S.iter) ' Wave #' num2str(iw) ' Probe # ' num2str(ip)])
                    meanProbeh = mean( abs(S.ProbeAmp{iw,ip}(S.bMask)).^2 );
                    ht = text(xlim(1), ylim(2), [' Mean Probe Int = ' num2str(meanProbeh, '%.1e')]);
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
                    set(gca,'xlim', xlim, 'ylim', ylim)
                    title(['Iteration # ' num2str(S.iter) ' Wave #' num2str(iw) ' Probe # ' num2str(ip) ' Amplitude'])
                    
                end
            end

        end % DisplayProbeAmp
        
        function [pl, hfig, hax] = DisplayOneProbeAmp(S, ip, varargin)
            if isempty(S.ProbeAmp),
                S.ReadReducedCube;
            end

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
        
        function [hfig, ha] = DisplayProbeCube(S, varargin)
            % [hfig, hax] = DisplayProbeCube(S, varargin)
            % ProbeCube is 2nd HDU of reduced data cube
            %
            %             CheckOption('hfig', [], varargin{:});
            %             CheckOption('iwv', ceil(S.NofW/2), varargin{:});
            %             CheckOption('blog', true, varargin{:});
            %             CheckOption('xlim', S.XlimDefault, varargin{:});
            %             CheckOption('ylim', S.YlimDefault, varargin{:});
            
            if isempty(S.ProbeModel),
                S.ReadProbeCube;
            end
            if isempty(S.bMask),
                S.ReadMaskCube;
            end
            
            hfig = CheckOption('hfig', [], varargin{:});
            iwvplot = CheckOption('iwv', ceil(S.NofW/2), varargin{:});
            bLog = CheckOption('blog', true, varargin{:});
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            
            [x, y] = CreateGrid(S.ProbeModel{1,1}, 1./S.ppl0);                        
            
            if bLog,
                funPlot = @(Ep) real(log10( abs(Ep).^2 ));
                %                 ImModel = real(log10(abs(S.ProbeModel{iwvplot,ip}).^2));
                %                 ImMeasure = real(log10(S.ProbeMeasAmp{iwvplot,ip}.^2));
                sctitle = 'log_{10} Norm Intensity';
            else
                funPlot = @(Ep) abs(Ep).^2;
                %                 ImModel = abs(S.ProbeModel{iwvplot,ip}).^2;
                %                 ImMeasure = S.ProbeMeasAmp{iwvplot,ip}.^2;
                sctitle = 'Linear Norm Intensity';
            end
                
            if isempty(hfig),
                hfig = figure_mxn(2,S.Nppair+1);
            else
                figure(hfig)
            end
            for ip = 1:S.Nppair,
                ha(ip) = subplot(2,S.Nppair+1,ip);
                imageschcit(x,y, funPlot(S.ProbeModel{iwvplot,ip})), axis image, colorbartitle(sctitle)
                set(gca,'xlim', xlim, 'ylim', ylim)
                xlabel('\lambda / D'), ylabel('\lambda / D')
                title(['wave #' num2str(iwvplot) ', it ' num2str(S.iter) ', Model Probe #' num2str(ip)])
                
                ha(S.Nppair+ip) = subplot(2,S.Nppair+1,S.Nppair+1+ip);
                imageschcit(x,y, funPlot(S.ProbeMeasAmp{iwvplot,ip})), axis image, colorbartitle(sctitle)
                set(gca,'xlim', xlim, 'ylim', ylim)
                xlabel('\lambda / D'), ylabel('\lambda / D')
                title(['wave #' num2str(iwvplot) ', it ' num2str(S.iter) ', Measure Probe #' num2str(ip)])
            end
            % match clim for all images
            climall = get(ha,'clim');
            if bLog, climmin = -9; else, climmin = 0; end
            clim = [climmin max([climall{:}])];
            set(ha,'clim',clim)

            
            % probe model is complex, take abs
            for ip = 1:S.Nppair,
                ProbeModelPlot{ip} = abs(S.ProbeModel{iwvplot,ip}).^2;
                ProbeMeasPlot{ip}  = abs(S.ProbeMeasAmp{iwvplot,ip}).^2;
                legstr{ip} = ['# ' num2str(ip)];
            end
            % plot limits for radial plot
            mean_rad = hypot(mean(xlim), mean(ylim));
            max_rad = max(hypot(xlim, ylim));
            radlim = mean_rad + (max_rad-mean_rad)*[-1 1];
            nr = min([128 ceil(2*(max_rad-mean_rad)*S.ppl0)]);
            %
            harad(1) = subplot(2,S.Nppair+1,S.Nppair+1);
            S.DisplayRadialPlot(ProbeModelPlot, 'dispradlim', radlim, 'nr', nr, 'hax', harad(1),'title', ['Iter #' num2str(S.iter) ', Model'], 'legstr', legstr);
            harad(2) = subplot(2,S.Nppair+1,2*(S.Nppair+1));
            S.DisplayRadialPlot(ProbeMeasPlot, 'dispradlim', radlim, 'nr', nr, 'hax',harad(2),'title', ['Iter #' num2str(S.iter) ', Measure'], 'legstr', legstr);
            
            ylim = get(harad,'ylim');
            set(harad,'ylim',[min([ylim{:}]), max([ylim{:}])])
            
        end % DisplayProbeCube
        
        function [hfig, hax] = DisplayCohInt(S, varargin)
            % [hfig, hax] = DisplayCohInt(S, varargin)
            % options:
            %             CheckOption('blog', true, varargin{:});
            %             CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            %             CheckOption('drawylimlines', [], varargin{:});
            %             CheckOption('clim', [], varargin{:});
            %             CheckOption('hax', [], varargin{:}); % put image on this axes
            %             CheckOption('xlim', [], varargin{:});
            %             CheckOption('ylim', [], varargin{:});
            
            if isempty(S.E_t),
                S.ReadReducedCube;
            end
            
            %             for iwv = 1:S.Nlamcorr,
            %                 CohInt{iwv} = abs(squeeze(S.E_t(iwv,:,:))).^2;
            %             end

            % options:
            % %  val = CheckOption(sOpt, valDefault, varargin)
            bLog = CheckOption('blog', true, varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            drawYlimLines = CheckOption('drawylimlines', [], varargin{:});
            clim = CheckOption('clim', [], varargin{:});
            haxuse = CheckOption('hax', [], varargin{:}); % put image on this axes
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            
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
            
            if isempty(clim),
                clim = pClim(pFun([S.CohInt{:}]));
            end
            
            [x, y] = CreateGrid(S.CohInt{1}, 1./(S.ppl0*S.PIAAMAG));

            for iwv = 1:S.Nlamcorr,
                if isempty(haxuse),
                    hax(iwv) = subplot(1, S.Nlamcorr, iwv);
                else
                    hax(iwv) =  haxuse(iwv);
                    axes(haxuse(iwv));
                end
              
                him(iwv) = imageschcit(x, y, pFun(S.CohInt{iwv})); axis image
                xlabel('\lambda / D')
                ylabel('\lambda / D')
                %caxis(clim);
                colorbartitle(cbartitle)                
                if ~isempty(S.NKTcenter), strlam = ['\lambda = ' num2str(S.NKTcenter(iwv)/S.NM) 'nm']; else strlam = ['Wave #' num2str(iwv)]; end
                title(['it#' num2str(S.iter) ', ' strlam])
            
            end
            set(hax,'xlim',xlim,'ylim',ylim,'clim',clim);
            
            DrawCircles(hax, drawRadii);
            DrawYlimLines(hax, drawYlimLines, drawRadii);
            
        end % DisplayCohInt
        
        function [hfig, haxlist, han] = DisplayAllInt(S, varargin)
            % display large table of unprobed, coh int, inc int images
            % uses:
            % S.DisplayImCubeUnProb
            % S.DisplayCohInt
            % S.DisplayIncInt
            % 
            % some options
            %    drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
            %    drawylimlines = CheckOption('drawylimlines', [], varargin{:})
            %    clim = CheckOption('clim', [-9 -6.5], varargin{:});

            bPlotRadialIntensity = CheckOption('DisplayRadialIntensity', true, varargin{:});
            hfig = CheckOption('hfig', [], varargin{:});
            
            % defaults that might be different
            varargin{end+1} = 'bLog'; varargin{end+1} = true;
            varargin{end+1} = 'clim'; varargin{end+1} = [-9 -6.5];
            
            % check that this instance is not empty
            if isempty(S.Nlamcorr)
                hfig = [];
                haxlist = [];
                return
            end %
            
            if isempty(hfig),
                hfig = figure;
            else
                % need to remove hfig from varargin
                iv = find(strcmp(varargin, 'hfig'));
                varargin{iv+1} = [];
            end
            
            %haxlist = zeros(3,S.Nlamcorr);
            if S.Nlamcorr == 1,
                nrow_ax = 1;
                ncol_ax = 3;
            else
                ncol_ax = S.Nlamcorr;
                nrow_ax = 3;
            end
            hfig = figure_mxn(hfig, nrow_ax, ncol_ax);
            
            % unprobed images
            figure(hfig);
            for ii = 1:S.Nlamcorr,
                haxlist(1,ii) = subplot(nrow_ax,ncol_ax,ii);
            end
            S.DisplayImCubeUnProb('hax',haxlist(1,:),varargin{:});
            % S.ImCubeUnProb{iwvpl}
            
            % Coh Int
            figure(hfig);            
            for ii = 1:S.Nlamcorr,
                haxlist(2,ii) = subplot(nrow_ax, ncol_ax,ii+S.Nlamcorr);
            end
            S.DisplayCohInt('hax',haxlist(2,:), varargin{:});
            % S.CohInt{iwv}
            
            % Inc Int
            figure(hfig);            
            for ii = 1:S.Nlamcorr,
                haxlist(3,ii) = subplot(nrow_ax, ncol_ax,ii+2*S.Nlamcorr);
            end
            S.DisplayIncInt('hax',haxlist(3,:), varargin{:});
            % S.IncIntEst{iwv}
            
            % make common clim
            %climlist = get(haxlist,'clim');
            %set(haxlist,'clim',[min([climlist{:}]) max([climlist{:}])])
            %set(haxlist,'clim',[-9 -6.5])

            % add text 'UnProbed', 'Modulated' and 'Unmodulated' to each row
            % Modulated:
            ylpos = get(get(haxlist(1,1),'YLabel'),'Position');
            han(1) = text(haxlist(1,1), ylpos(1) - 2, ylpos(2), 'UnProbed' ...
                , 'Rotation',90 ...
                ,'HorizontalAlignment','center' ...
                ,'VerticalAlignment','top' ...
                ,'VerticalAlignment','bottom' ...
                ... ,'Position', ylpos - [2 0 0] ...
                ,'FontSize', 18 ...
                ,'Color','b' ...
                ,'FontWeight','bold' ...
                );
            ylpos = get(get(haxlist(2,1),'YLabel'),'Position');
            han(2) = text(haxlist(2,1), ylpos(1) - 2, ylpos(2), 'Modulated' ...
                , 'Rotation',90 ...
                ,'HorizontalAlignment','center' ...
                ,'VerticalAlignment','top' ...
                ,'VerticalAlignment','bottom' ...
                ... ,'Position', ylpos - [2 0 0] ...
                ,'FontSize', 18 ...
                ,'Color','b' ...
                ,'FontWeight','bold' ...
                );
            ylpos = get(get(haxlist(3,1),'YLabel'),'Position');
            han(3) = text(haxlist(3,1), ylpos(1) - 2, ylpos(2), 'UnModulated' ...
                , 'Rotation',90 ...
                ,'HorizontalAlignment','center' ...
                ,'VerticalAlignment','top' ...
                ,'VerticalAlignment','bottom' ...
                ... ,'Position', ylpos - [2 0 0] ...
                ,'FontSize', 18 ...
                ,'Color','b' ...
                ,'FontWeight','bold' ...
                );

            if bPlotRadialIntensity,
                % for radial plots, need to set bMask to exactly the
                % control region
                
                S.DisplayRadialIntensity(varargin{:});
            end
            
        end
        
        function [hfig, hax, hl, rplot, IntRad] = DisplayRadialIntensity(S, varargin)
            % [hfig, hax, hl, rplot, IntRad] = DisplayRadialIntensity(S, varargin)
            %
            % radial plot of full band mean total, unmodulated, modulated
            
            xlim = CheckOption('xlim', [], varargin{:});
            ylim = CheckOption('clim', [], varargin{:}); % same as clim for the images
            ylim = 10.^ylim; % assumes semilogy
            %clim = CheckOption('clim', [-9 -6.5], varargin{:});
            drawylimlines = CheckOption('drawylimlines', [], varargin{:});
            
            if isempty(S.ImCubeUnProbFullBand),
                S.ReadImageCube;
            end
            if isempty(S.CohIntFullBand) || isempty(S.IncIntEstFullBand),
                S.ReadReducedCube;
            end
            
            if ~isempty(drawylimlines)
                [x, y, X, Y, R] = CreateGrid(S.ImCubeUnProbFullBand, 1./(S.ppl0*S.PIAAMAG));
                bMaskuse = S.bMask & Y > drawylimlines(1); % & (~squeeze(S.bMask_badpix(1,:,:)));
            else
                bMaskuse = S.bMask;
            end
            [sHandles, rplot, IntRad] = S.DisplayRadialPlot( ...
                {S.ImCubeUnProbFullBand, S.CohIntFullBand, S.IncIntFullBand}, ...
                'legstr', {'Total','Modulated','Unmodulated'}, ...
                'plotmean', false, 'bMask', bMaskuse, varargin{:}, 'nr', 88);
            set(sHandles.hline,'LineWidth',2)
            set(sHandles.hline(1),'color','k')
            set(sHandles.hline(3),'color','b')
            if ~isempty(xlim), set(sHandles.hax,'xlim',xlim), end
            if isempty(ylim)
                ylim = get(sHandles.hax,'ylim'); set(sHandles.hax,'ylim', [max([ylim(1) 1e-9]), ylim(2)]);
            else
                set(sHandles.hax,'ylim',ylim);
            end
            
    
            
        end % DisplayRadialIntensity
        
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
            %clim = CheckOption('clim', [], varargin{:});
            
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
        
        function [hfig, hax] = DisplayEfields(S, varargin)
            % [hfig, ha] = DisplayEfields(S)
            % 
            % testbed measured E fields
            % amp, real, imag x each wavelength
            %
            %             xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            %             ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            %   CheckOption('hfig', [], varargin{:});
            %   CheckOption('clim', [], varargin{:});
            %   CheckOption('drawradii', S.DrawradiiDefault, varargin{:});
                        
            if isempty(S.E_t),
                S.ReadReducedCube;
            end
            
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            hfig = CheckOption('hfig', [], varargin{:});
            clim = CheckOption('clim', [], varargin{:});
            drawRadii = CheckOption('drawradii', S.DrawradiiDefault, varargin{:});

            % title
            sRI = ['run #' num2str(S.runnum) ', iter #' num2str(S.iter)];
            
            %
            [nw, nr, nc] = size(S.E_t);
            if nw ~= S.NofW, error('size E_t does not match'); end
            [x, y] = CreateGrid([nc nr], 1./S.ppl0);

            %
            Nplr = 3; % # of rows
            if isa(hfig,'matlab.ui.Figure'), figure(hfig)
            else, hfig = figure_mxn(Nplr,S.NofW);
            end
            hax = zeros(Nplr, S.NofW);
            
            % first row is amplitude
            for iwv = 1:S.NofW,
                hax(1,iwv) = subplot(Nplr,S.NofW,iwv);
                imageschcit(x, y, abs(squeeze(S.E_t(iwv,:,:)))), colorbar
                title(['Amp, Wave#' num2str(iwv) ', ' sRI])
            end
            
            % second row is real part
            for iwv = 1:S.NofW,
                hax(2,iwv) = subplot(Nplr, S.NofW, S.NofW + iwv);
                imageschcit(x, y, real(squeeze(S.E_t(iwv,:,:)))), colorbar
                title(['Real, Wave#' num2str(iwv) ', ' sRI])
            end
            
            % third row is imag part
            for iwv = 1:S.NofW,
                hax(3,iwv) = subplot(Nplr, S.NofW, 2*S.NofW + iwv);
                imageschcit(x, y, imag(squeeze(S.E_t(iwv,:,:)))), colorbar
                title(['Imag, Wave#' num2str(iwv) ', ' sRI])
            end

            %
            set(hax,'xlim',xlim,'ylim',ylim)
            
            % clim
            if isempty(clim),
                Euse = S.E_t(abs(S.E_t(:)) > 0);
                aclim = AutoClim(abs(Euse), 'one-sided', true);
                clim  = AutoClim([real(Euse) imag(Euse)], 'symmetric', true);
            else
                aclim = [0 clim(2)];                
            end
            set(hax(1,:),'clim',aclim)  % abs plots
            set(hax(2:3,:),'clim',clim) % real, imag plots
  
            DrawCircles(hax, drawRadii);
            
        end % DisplayEfields

        function [hfig, ha, sMetrics] = DisplayDEfields(S, Sref, varargin)
            % [hfig, ha, sMetrics] = DisplayDEfields(S, Sref, varargin)
            %
            % 4 x NofW, dE_t real, imag, dE_m real, imag
            %
            % sMetrics = 
            %                     'type', 'dEfields' ...
            %                     ,'rmsdE_t', nan ...
            %                     ,'rmsdE_m', nan ...
            %                     ,'dE_t', [] ...
            %                     ,'dE_m', [] ...
            % 
            % CheckOption('hfig', [], varargin{:});
            % CheckOption('clim', [], varargin{:});
            % CheckOption('nodisplay', false, varargin{:});
            % CheckOption('xlim', S.XlimDefault, varargin{:});
            % CheckOption('ylim', S.YlimDefault, varargin{:});

            hfig = CheckOption('hfig', [], varargin{:});
            clim = CheckOption('clim', [], varargin{:});
            bNodisplay = CheckOption('nodisplay', false, varargin{:});
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            
            if isempty(S.E_t),
                S.ReadReducedCube;
            end

            if isempty(Sref.E_t),
                Sref.ReadReducedCube;
            end

            [nw, nr, nc] = size(S.E_t);
            if nw ~= S.NofW, error(['number of wavelengths inconsistent']); end
            
            % check S and Sref match
            % check that S.E_t and Sref.E_t are same size
            % also catches if one is no data
            if ~isequal(size(S.E_t), size(Sref.E_t)),
                disp(['iter ' num2str(S.iter) ' and iter ' num2str(Sref.iter) ' mismatch, skipping DisplayCEfields']);
                [hfig, ha] = deal([]);

                sMetrics = struct(...
                    'type', 'dEfields' ...
                    ,'rmsdE_t', nan ...
                    ,'rmsdE_m', nan ...
                    ,'dE_t', [] ...
                    ,'dE_m', [] ...
                );

                return
            end

            %sRI = ['run #' num2str(S.runnum) ', iter #' num2str(S.iter) '--' num2str(Sref.iter)];
            sRI = ['\DeltaE Iter #' num2str(S.iter) ' - ' num2str(Sref.iter)];
            %sRI = '';
                        
            %
            dE_t = S.E_t - Sref.E_t;
            dE_m = S.E_m - Sref.E_m;
            
            % calculate dE metrics
            % pixels to use for metrics
            bMaskuse = S.bMaskSc & Sref.bMaskSc;
            [rmsdE_t, rmsdE_m] = deal(zeros(S.NofW,1));
            for iwv = 1:S.NofW,
                % use only score region for rms dE
                bMaskiwl = bMaskuse & (S.IncInt{iwv} >= 0) & (Sref.IncInt{iwv} >= 0);
                rmsdE_t(iwv) = sqrt(mean(abs(dE_t(iwv,bMaskiwl)).^2));
                rmsdE_m(iwv) = sqrt(mean(abs(dE_m(iwv,bMaskiwl)).^2));
            end
            sMetrics = struct(...
                'type', 'dEfields' ...
                ,'rmsdE_t', rmsdE_t ...
                ,'rmsdE_m', rmsdE_m ...
                ,'dE_t', dE_t ...
                ,'dE_m', dE_m ...            
                );
            
            % if no display, return metrics and skip graphs
            if bNodisplay,
                hfig = [];
                ha = [];
                return
            end
            
            
            % top row = real(DE_t)
            % 2nd row = imag(DE_t)
            % 3rd row = real(DE_m)
            % 4ty row = imag(DE_m)

            % prepare figure
            Nplr = 4;
            if isa(hfig,'matlab.ui.Figure'),
                figure(hfig)
            else,
                hfig = figure_mxn(Nplr,S.NofW);
            end
            [x, y] = CreateGrid([nr nc], 1./S.ppl0);
            ha = zeros(Nplr,S.NofW);
            for iwv = 1:S.NofW,
                % subplot #
                iptr = iwv+0*S.NofW;
                ipti = iwv+1*S.NofW;
                ipmr = iwv+2*S.NofW;
                ipmi = iwv+3*S.NofW;
                
                % change in unprobed image intensity
                ha(1,iwv) = subplot(Nplr, S.NofW, iptr);
                imageschcit(x,y,squeeze(real(dE_t(iwv,:,:)))); %colorbar
                title(['Measure: real{\DeltaE}, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])
                
                ha(2,iwv) = subplot(Nplr, S.NofW, ipti);
                imageschcit(x,y,squeeze(imag(dE_t(iwv,:,:)))); %colorbar
                title(['Measure: imag{\DeltaE}, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])

                ha(3,iwv) = subplot(Nplr, S.NofW, ipmr);
                imageschcit(x,y,squeeze(real(dE_m(iwv,:,:)))); %colorbar
                title(['Model: real{\DeltaE}, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])
                
                ha(4,iwv) = subplot(Nplr, S.NofW, ipmi);
                imageschcit(x,y,squeeze(imag(dE_m(iwv,:,:)))); %colorbar
                title(['Model: imag{\DeltaE}, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])
                
            end            

            % xlim, ylim
            set(ha,'xlim',xlim,'ylim',ylim)

            % clim
            if isempty(clim),
                dEtmp = [dE_t dE_m];
                dEuse = dEtmp(abs(dEtmp(:)) > 0);
                clim = AutoClim([real(dEuse) imag(dEuse)],'symmetric',true);
                %set(ha,'clim',median(climE))
            else
                %set(ha,'clim',clim)
            end
            set(ha,'clim',clim)

            % put colorbars on the right-most axes
            for ii = 1:4,
                colorbar('peer',ha(ii,end))
            end
            
            % make a title at the top
            % delete previous annotation if it exists
            hantmp = findobj(gcf,'type','textboxshape');
            if ~isempty(hantmp)
               hantmp.delete; 
            end
            % new annotation
            han = annotation('textbox', [0.5 0.8 0.2 0.2], 'String', sRI, ...
                'FitBoxToText', 'on', 'LineStyle', 'none', ...
                'FontSize', 24, 'Color', 'r', 'FontWeight', 'bold');
            set(han,'HorizontalAlignment','center')
            % center horizontally
            ppp = get(han,'Position');
            set(han,'Position',[0.5 - 0.5*ppp(3) ppp(2:end)])
            % so it can be found and deleted later
            set(get(han,'parent'),'HandleVisibility','on')

        end % DisplayDEfields

        function [hfig, ha, sCmetrics] = DisplayCEfields(S, Sref, varargin)
            % [hfig, ha] = DisplayCEfields(S, Sref, varargin)
            % correlation metrics DE_t .* conj(DE_m)
            %
            %             CheckOption('hfig', [], varargin{:});
            %             CheckOption('clim', [], varargin{:});
            %             CheckOption('PSF_thresh_nsig', 4, varargin{:});
            %             CheckOption('debug', false, varargin{:});
            %             CheckOption('bMaskDisplay', [], varargin{:}); % default is mask from CohInt
            %             CheckOption('nodisplay', false, varargin{:}); % calc metrics and return, don't display graph
            %             xlim = CheckOption('xlim', S.XlimDefault, varargin{:});
            %             ylim = CheckOption('ylim', S.YlimDefault, varargin{:});

            hfig = CheckOption('hfig', [], varargin{:});
            clim = CheckOption('clim', [], varargin{:});
            PSF_thresh_nsig = CheckOption('PSF_thresh_nsig', 4, varargin{:});
            bDebugAutoMetric = CheckOption('debug', false, varargin{:});
            bMaskDisplay = CheckOption('bMaskDisplay', [], varargin{:});
            bNodisplay = CheckOption('nodisplay', false, varargin{:});
            xlim = CheckOption('xlim', S.XlimDefault, varargin{:}); 
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});
            
            if isempty(S.E_t),
                S.ReadReducedCube;
            end

            % validate Sref
            if ~isa(Sref, 'CRunData'),
                Sref = CRunData(S.runnum, Sref);
            end
            
            if isempty(Sref.E_t),
                Sref.ReadReducedCube;
            end

            [nw, nr, nc] = size(S.E_t);
            if nw ~= S.NofW, error(['number of wavelengths inconsistent']); end

            % check that S.E_t and Sref.E_t are same size
            % also catches if one is no data
            if ~isequal(size(S.E_t), size(Sref.E_t))
                disp(['iter ' num2str(S.iter) ' and iter ' num2str(Sref.iter) ' mismatch, skipping DisplayCEfields']);
                [hfig, ha] = deal([]);

                sCmetrics = struct(...
                    'type', 'CEfields' ...
                    ,'CP', nan ...
                    ,'CC', nan ...
                    ,'angle_CC', nan ...
                    ,'mag_dEm_dEt', nan ...
                    ,'mse', nan ...
                    ,'CP_definition', ' <dEt.dEm>/<dEm.dEm> ' ...
                    ,'CC_definition', ' <dEt.dEm>/sqrt(<dEt.dEt><dEm.dEm>) ' ...
                    ,'mag_dEm_dEt_definition', ' sqrt(<dEm.dEm>/<dEt.dEt>) ' ...
                    ,'mse_definition', ' mean(abs(dEt - dEm).^2) ' ...
                );

                return
            end
            
            % title string
            sRI = ['iter #' num2str(S.iter) '--' num2str(Sref.iter)];
            
            % calculate correlation metrics
            dE_t = S.E_t - Sref.E_t;
            dE_m = S.E_m - Sref.E_m;
            CE   = conj(dE_m) .* dE_t ./sqrt( (dE_m(:)'*dE_m(:)).*(dE_t(:)'*dE_t(:)) );
            
            % from Joon:
            % CP =  <dEt.dEm>/<dEm.dEm>
            % CC =  <dEt.dEm>/sqrt(<dEt.dEt><dEm.dEm>)
            % |DE| = sqrt(<dEm.dEm>/<dEt.dEt>)
            bMaskCube = false(size(S.E_t));
            [nwtmp, nrtmp, nctmp] = size(S.E_t);
            if nwtmp ~= S.NofW, error('cube size wrong'); end
            for ii = 1:S.NofW,
                bMaskCube(ii,:,:) = S.bMaskSc & (S.IncInt{ii} > 0);
            end
            dE_mu = dE_m(bMaskCube);
            dE_tu = dE_t(bMaskCube);

            sCmetrics = struct(...
                'type', 'CEfields' ...
                ,'CP', (dE_mu'*dE_tu)./(dE_mu'*dE_mu) ...
                ,'CC', (dE_mu'*dE_tu)./sqrt( (dE_mu'*dE_mu).*(dE_tu'*dE_tu) ) ...
                ,'angle_CC', angle((dE_mu'*dE_tu)) ...
                ,'mag_dEm_dEt', sqrt( (dE_mu'*dE_mu)./(dE_tu'*dE_tu) ) ...
                ,'mse', mean( abs(dE_tu - dE_mu).^2 ) ...
                ,'CP_definition', ' <dEt.dEm>/<dEm.dEm> ' ...
                ,'CC_definition', ' <dEt.dEm>/sqrt(<dEt.dEt><dEm.dEm>) ' ...
                ,'mag_dEm_dEt_definition', ' sqrt(<dEm.dEm>/<dEt.dEt>) ' ...
                ,'mse_definition', ' mean(abs(dEt - dEm).^2) ' ...
                );
            
            if bNodisplay,
                hfig = [];
                ha = [];
                return
            end
            
            % prepare display figure
            Nplr = 2;
            if isa(hfig,'matlab.ui.Figure'),
                figure(hfig)
            else
                hfig = figure_mxn(Nplr,S.NofW);
            end
            [x, y] = CreateGrid([nr nc], 1./S.ppl0);
            ha = zeros(Nplr,S.NofW);

            for iwv = 1:S.NofW,
                
                % correlation amplitude for this iwv
                CEampiwv = squeeze(abs(CE(iwv,:,:)));

                % to make the phase plot cleaner, only plot phase where
                % there are speckles. 
                if isempty(bMaskDisplay),
                    % Let's make a mask from Coh Intensity (show where
                    % there are bright speckles)
                    [~, bMaskce] = AutoMetric(S.CohInt{iwv}, [], ...
                        struct('image_type','psf','logPSF', true, ...
                        'debug',bDebugAutoMetric,'PSF_thresh_nsig',PSF_thresh_nsig));
                    
                    %                     % Or, Let's make a mask from speckle
                    %                     % correlation intensity
                    %                     [~, bMaskce] = AutoMetric(CEampiwv, [], ...
                    %                         struct('image_type','psf','logPSF',true,...
                    %                         'debug',bDebugAutoMetric,'PSF_thresh_nsig',PSF_thresh_nsig));
                else
                    bMaskce = bMaskDisplay;
                end
                
                % add to display mask pixels where amp == 0 so pha is not
                % displayed
                bMaskce(CEampiwv <= 0) = false;
                
                % apply display mask
                CEphaiwv = squeeze(angle(CE(iwv,:,:)));
                CEphaiwv(~bMaskce) = nan;
                
                % plot log amp
                CEampiwvlog = real(log10(CEampiwv));
                CEampiwvlog(~bMaskce) = nan;
                
                % subplot #
                iamppl = iwv+0*S.NofW;
                iphapl = iwv+1*S.NofW;
                
                % amplitude of projection
                
                ha(1,iwv) = subplot(Nplr, S.NofW, iamppl);
                imageschcit(x,y, CEampiwvlog); colorbar
                title([sRI ', CE = |dE_m''dE_t|/\surd{<dE_t.dE_t><dE_m.dE_m>}, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])
                
                ha(2,iwv) = subplot(Nplr, S.NofW, iphapl);
                him_pha = imageschcit(x,y,CEphaiwv/pi);                 
                %mesh(x, y, CEphaiwv/pi, 'linestyle','none'); view(2); % doesn't plot NaN
                title([sRI ', \angle{CE}, ' num2str(S.NKTcenter(iwv)/S.NM) 'nm'])
                colorbartitle('Phase (\pi rad)')
                % circular colormap for phase plots
                set(ha(2,iwv),'clim',[-1 1]) % pi radians
                colormap(ha(2,iwv), hsv) %phasemap) % hsv also works
                set(him_pha, 'AlphaData', ~isnan(CEphaiwv))
                %axis(ha(2,iwv), 'image')
            end

            % xlim, ylim
            set(ha,'xlim',xlim,'ylim',ylim)

            % make clim for abs plots the same for all wavelengths
            if isempty(clim)
                if S.NofW > 1,
                    cclim = get(ha(1,:),'clim'); % returns a cell array
                    clim = [min([cclim{:}]) max([cclim{:}])];
                else
                    clim = get(ha(1,1),'clim');
                end
            end
            set(ha(1,:),'clim',clim)
            
        end % DisplayCEfields

        function [hfig, ha] = DisplayPampzero(S, varargin)
            % [hfig, ha] = DisplayPampzero(S, varargin)

            xlim = CheckOption('xlim', S.XlimDefault, varargin{:}); 
            ylim = CheckOption('ylim', S.YlimDefault, varargin{:});

            [nw, ny, nx] = size(S.bPampzero);
            [x, y] = CreateGrid([nx ny], 1./S.ppl0);

            hfig = figure_mxn(1,S.NofW);
            for iw = 1:S.NofW,
               ha(iw) = subplot(1,S.NofW,iw);
               
               imageschcit(x, y, squeeze(S.bPampzero(iw,:,:))), axis image
               xlabel('\lambda / D')
               ylabel('\lambda / D')
               
               set(gca,'xlim', xlim,'ylim', ylim)
               
            end
            
        end % DisplayPampzero

        function [hfig, hax, sMetrics] = DisplayDMv(S, dmvref, varargin)
            % [hfig, hax, sMetrics] = S.DisplayDMv([], varargin)
            % [hfig, hax, sMetrics] = S.DisplayDMv(Sref, varargin)
            % [hfig, hax, sMetrics] = S.DisplayDMv({refDM1v_fits, refDM2v_fits}, varargin)
            %
            % CheckOption('climdelta', [], varargin{:});
            % CheckOption('hfig', [], varargin{:});
            
            if nargin < 2, dmvref = []; end
            
            if isempty(S.DMvCube)
                S.ReadDMvCube;
            end
            
            % initialize return values
            [hfig, hax] = deal([]);
            sMetrics = struct(...
                'type', 'DMv' ...
                , 'rmsdDMv', NaN ...
                , 'dDMv', {} ...
                );

            % check if empty instance
            if isempty(S.DMvCube),
                return
            end
            
            climDelta = CheckOption('climdelta', [], varargin{:});
            hfig = CheckOption('hfig', [], varargin{:});

            % extract the DV v to plot
            for idm = 1:S.Ndm,
                DMvtmp = squeeze(S.DMvCube{idm}(:,:,1));
                rmsDMv(idm) = rms(DMvtmp(DMvtmp>0));
                
                DMv{idm} = DMvtmp;
                strDM{idm} = ['it#' num2str(S.iter) ...
                    ', DM' num2str(idm) ...
                    ', ' num2str(rmsDMv(idm),'%.4f') 'V rms'];
            end

            % is there a reference DMv
            refDMv = [];
            strRefDM = [];

            if ~isempty(dmvref),
                if isa(dmvref, 'cell')
                    [refDM1v_fn, refDM2v_fn] = deal(dmvref{:});
                    
                    refDMv{1} = fitsread(PathTranslator(refDM1v_fn));
                    aatmp = regexp(refDM1v_fn, '/', 'split');
                    strRefDM{1} = pwd2titlestr(aatmp{end});
                                                    
                    refDMv{2} = fitsread(PathTranslator(refDM2v_fn));
                    aatmp = regexp(refDM2v_fn, '/', 'split');
                    strRefDM{2} = pwd2titlestr(aatmp{end});
                elseif isa(dmvref, 'CRunData')
                    Sref = dmvref;
                    if isempty(Sref.DMvCube)
                        Sref.ReadDMvCube;
                    end
                    
                    for idm = 1:Sref.Ndm,
                        refDMv{idm} = squeeze(Sref.DMvCube{idm}(:,:,1));
                        strRefDM{idm} = ['it#' num2str(Sref.iter) ', DM' num2str(idm)];
                    end
                    
                end % if isa(dmvref
                
            end % if ~isempty(dmvref)
                    
            Nr = 1 + ~isempty(refDMv);
            if isempty(hfig),
                hfig = figure_mxn(Nr,S.Ndm);
            else
                figure(hfig);
            end
            
            rmsdDMv = zeros(1,S.Ndm);
            for idm = 1:S.Ndm,

                % plot S DMv
                hax(idm) = subplot(Nr, S.Ndm, idm);
                imageschcit(0,0,DMv{idm}), 
                colorbartitle('Vmu')
                title(strDM{idm})
                
                % if Ref defined:
                if ~isempty(refDMv),
                    hax(idm+S.Ndm) = subplot(Nr, S.Ndm, idm+S.Ndm);
                
                    dDMv = DMv{idm} - refDMv{idm};
                    rmsdDMv(idm) = rms(dDMv(abs(dDMv)>0));

                    imageschcit(0,0,dDMv)
                    colorbartitle('Vmu')
                    title(['\Delta ' strRefDM{idm} ', ' num2str(rmsdDMv(idm),'%.4f') 'V rms'])
                    
                    % save
                    cdDMv{idm} = dDMv;
                end
                
            end % for idm
            
            % equalize clim for DMv
            if ~isempty(refDMv),
                aclim = AutoClim([DMv{:}],'one-sided',true);
                set(hax(1:S.Ndm),'clim',aclim);
                
            end % refDMv

            % equalize clim for ddm
            if ~isempty(refDMv),
                if ~isempty(climDelta),
                    set(hax(S.Ndm+1:end),'clim',climDelta)
                elseif range(reshape([cdDMv{:}],[],1)) > 0,
                    aclim = AutoClim([cdDMv{:}],'symmetric',true);
                    set(hax(S.Ndm+1:end),'clim',aclim);
                end                
            end % refDMv
            
            sMetrics = struct(...
                'type', 'DMv' ...
                ,'rmsdDMv', rmsdDMv ...
                ,'dDMv', {cdDMv} ... % {} so that the cell array is assigned to one struct
                );
            
            %             fprintf('rms dDMv1 = %.3f Vmu\n',rmsdDMv1);
            %             fprintf('rms dDMv2 = %.3f Vmu\n',rmsdDMv2);
            %
        end % DisplayDMv
        
        function [hfig, hax] = DisplayDMvProbe(S, varargin)
            % [hfig, hax] = DisplayDMvProbe(S, varargin)
            
            if isempty(S.DMvCube),
                S.ReadDMvCube;
            end
            
            idm = CheckOption('idm', 1, varargin{:}); % which DM is probing?
            
            % probed DM:
            DMv = S.DMvCube{idm};
            
            [nacty, nactx, nsli] = size(DMv);
            npr = nsli-1;
            
            hfig = figure_mxn(2,npr/2);
            for ii = 1:npr/2,
                isl = 2*ii; % slice into dmv cube
                hax(ii) = subplot(2,npr/2,ii);
                imageschcit(0,0,squeeze(DMv(:,:,isl)-DMv(:,:,1)))
                title(['iter #' num2str(S.iter) '; Probe #' num2str(ii) '; +ve'])
                colorbartitle('Vmu')
                
                hax(ii+npr/2) = subplot(2,npr/2,ii+npr/2);
                imageschcit(0,0,squeeze(DMv(:,:,isl+1)-DMv(:,:,1)))
                title(['iter #' num2str(S.iter) '; Probe #' num2str(ii) '; -ve'])
                colorbartitle('Vmu')
                                    
            end % for each pair
            
            % make uniform clim
            clim = get(hax,'clim'); % a cell array
            set(hax,'clim', max(abs([max([clim{:}]) min([clim{:}])]))*[-1 1])
            
        end % DisplayDMvProbes
        
        function [hfig, hax] = DisplayImCubeDelProbes(S, iwv, varargin)
            % [hfig, hax] = DisplayImCubeDelProbes(S, iwv, varargin)
            % 
            % S.ImCubeDelProb{iwv,ipr} = ImPrPlus - ImPrMinus;
            %
            % options:
            
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

        function Sppt = GenerateReportPpt(S, varargin)
            % Sppt = GenerateReportPpt(S, varargin)
            %
            %      Sppt = CheckOption('Sppt', S.Sppt, varargin{:});
            %      ppt_fn = CheckOption('fn', '', varargin{:});
            %      Sref = CheckOption('ref', [], varargin{:});
        
            Sppt = CheckOption('Sppt', S.Sppt, varargin{:});
            ppt_fn = CheckOption('fn', '', varargin{:});
            Sref = CheckOption('ref', [], varargin{:});
            
            if isempty(Sppt)
                Sppt = Cppt(ppt_fn);
            end
            
            function CopyFigSlide(hfig)
                set(hfig,'Position',0.6*get(hfig,'Position'));
                Sppt.CopyFigNewSlide(hfig);
            end
            
            hfig = S.DisplayAllInt;
            CopyFigSlide(hfig);
            
            hfig = S.DisplayProbeAmp;
            CopyFigSlide(hfig);
            
            hfig = S.DisplayEfields;
            CopyFigSlide(hfig);
                        
            % Displays that require a reference
            if ~isempty(Sref)
                if isnumeric(Sref),
                    Sref = CRunData(S.runnum, Sref);
                end
                
                hfig = S.DisplayDEfields(Sref);
                CopyFigSlide(hfig);
                
                hfig = S.DisplayCEfields(Sref);
                CopyFigSlide(hfig);
                
                hfig = S.DisplayDMv(Sref);
                CopyFigSlide(hfig);
                
            end % if reference
            
            S.Sppt = Sppt;
            
        end % GenerateReportPpt
        
        
        
    end % methods
        
end % classdef

% utilities
function val = CheckOption(varstring, defaultval, varargin)
% val = CheckOption(varstring, defaultval, varargin)
%
% utility for checking varargin for an option
% varstring = string keyword
% if varstring is found in varargin{:},
% the value of the following entry in varargin is returned
% else the defaultval is returned
%
% example:
% RminSc = CheckOption('RminSc', S.RminSc, varargin{:});

iv = find(strcmpi(varargin, varstring));

if ~isempty(iv),
    val = varargin{iv+1};
else
    val = defaultval;
end

end % CheckOption

function DrawCircles(hax, drawRadii)
% DrawCircles(hax, drawRadii)

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
% DrawThetaLines(hax, drawTheta, drawRadii)

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

function DrawYlimLines(hax, drawYlimLines, drawRadii)
% DrawYlimLines(hax, drawYlimLines, drawRadii)

if isempty(drawYlimLines),
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
    for ith = 1:length(drawYlimLines),
        yline = drawYlimLines(ith);
        x0 = sqrt(r1.^2 - yline.^2);
        if abs(yline) > r0
            hl = plot([x0 -x0], yline*[1 1], '-r');
        else
            x1 = sqrt(r0.^2 - yline.^2);
            hl = plot([-x0 -x1], yline*[1 1], '-r', [x1 x0], yline*[1 1], '-r');
        end
        set(hl,'linewidth',1)
        
    end % for each y line
    hold off;

end % for each axes

end % DrawThetaLines
