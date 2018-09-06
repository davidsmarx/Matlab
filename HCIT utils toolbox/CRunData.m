classdef CRunData < handle & CConstants
% class CRunData
%
% S = CRunData(runnum, iter, sOptions)
%     default  sOptions = struct(...
%                 'results_pn', ['../results/run' num2str(runnum, '%03d') '/'] ...
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
% [hfig, hax] = DisplayImCubeUnProb(S, sOptions)
%             sOpt = struct(...
%                 'bPlotLog', false ...
%                 ,'drawRadii', [] ...
%                 ,'clim', [] ...
%                 );
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
%    'drawradii'
%    'bLog'
%
% [hfig, hax] = DisplayDMv(S)
% [hfig, hax] = DisplayDMv(S, Sref)
% [hfig, hax] = DisplayDMv(S, DM1ref_fits_fn, DM2ref_fits_fn)

%    'bLog', true or false (default = false)
%    'drawradii', [list of radii]
%
% [hfig, ha] = DisplayDEfields(S, Sref, hfig)

    properties
        
        % SPC disc:
        lambda          = [542 553 565 577 588]*CConstants.NM;
        %ilamcorr        = [0, 1, 2, 3, 4];
        ppl0            = 4.01;
        RminSc          = 6.3; % lam/D
        RmaxSc          = 19.5;
        Nbscan          = 6;

        Reduced_fn = '';
        Rundir_fn  = '';
        Results_pn = '';
        
        runnum
        iter
        
        % results:
        IncInt
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
        ImCubeUnProb
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
        
        DMvCube
        ProbeModel
        ProbeMeasAmp
        ProbeMeasCross
        ProbeRes
        ProbeCon
        
        rplot
        IntRad
        
    end % properties
        
    % reduced fits:
    % Primary: reducedcube
    %   Image: maskcube
    %   Image: probecube
    
    methods
        function S = CRunData(runnum, iter, sOptin)
            
            if ~exist('sOptin','var'), sOptin = struct(); end
            % validate sOptions, default values:
            sOptions = struct(...                  %'results_pn', ['../results/run' num2str(runnum, '%03d') '/'] ...
                'results_pn', PathTranslator(['/home/dmarx/ln_dst_data/EFC/HLC/run' num2str(runnum, '%03d') '/']) ...
                ,'PCtemp_pn', 'C:\Users\dmarx\HCITdatatemp\' ...
                ,'debug',false ...
                );
            % replace default values with given options
            ffields = fieldnames(sOptin);
            for ii = 1:length(ffields),
                sOptions.(ffields{ii}) = sOptin.(ffields{ii});
            end

            % default if not already defined by options:
            if ~isfield(sOptions,'Rundir_pn')
                sOptions.Rundir_pn =  [sOptions.results_pn 'rundir/'];
            end
            if ~isfield(sOptions,'Reduced_pn')
                sOptions.Reduced_pn = [sOptions.results_pn 'rundir/reduced/'];
            end

            % if only relative paths are used, don't need PathTranslator
            %results_pn = ['../results/run' num2str(runnum, '%03d') '/'];
            if ~exist(PathTranslator(sOptions.results_pn), 'dir'),
                error(['not found: ' sOptions.results_pn]);
            end
            S.Results_pn = sOptions.results_pn;

            % array of iteration #'s
            if ~exist('iter','var') || isempty(iter),
                error('must define iterv iterations');
            end
     
            % unzip the fits files if necessary
            % basename, such as 'run603it00000.fits'
            s_bn = ['run' num2str(runnum,'%03d') 'it' num2str(iter,'%05d') '.fits'];
            % 'local' fits files are on PC, if ispc
            if ispc,                
                if ~exist([sOptions.PCtemp_pn 'reduced\'],'dir')
                    mkdir([sOptions.PCtemp_pn 'reduced\']);
                end
                sReduced_fn = [sOptions.PCtemp_pn 'reduced\' s_bn];
            else
                sReduced_fn = [sOptions.Reduced_pn s_bn];
            end
            
            % first check if 'local' fits files exist, then skip ahead
            if ~exist(PathTranslator(sReduced_fn),'file'),
                % source location:
                sReduced_gz = [sOptions.Reduced_pn s_bn '.gz'];
                if ~exist(PathTranslator(sReduced_gz), 'file'),
                    warning(['data file: ' sReduced_fn '.gz' ' does not exist, no reduced data available']);
                    sReduced_fn = [];
                else
                    
                    % if desktop pc, faster to copy local then unzip
                    if ispc,
                        copyfile(PathTranslator(sReduced_gz), [sReduced_fn '.gz']);
                        sReduced_gz = [sReduced_fn '.gz'];
                    end
                end
                % sReduced_gz is now either on server or local pc
                gunzip(sReduced_gz);
            end
            
            % repeat local logic for rundir
            if ispc,
                sRundir_fn = [sOptions.PCtemp_pn s_bn];
            else,
                sRundir_fn = PathTranslator([sOptions.Rundir_pn s_bn]);
            end
            
            % if 'local' rundir fits file does not exist, unzip
            if ~exist(sRundir_fn,'file'),
                % source location:
                sRundir_gz = [sOptions.Rundir_pn s_bn '.gz'];
                if ~exist(PathTranslator(sRundir_gz),'file'),
                    error(['data file: ' sRundir_gz ' does not exist']);
                end
                
                if ispc,
                    copyfile(PathTranslator(sRundir_gz), [sRundir_fn '.gz']);
                    sRundir_gz = [sRundir_fn '.gz'];
                end
                    
                % sRundir_gz is now either local pc or still server
                % gunzip should create sRundir_fn
                gunzip(sRundir_gz);
                
            end
            
            finfo = fitsinfo(sRundir_fn);
            S.ImKeys = finfo.PrimaryData.Keywords;
            
            if ~isempty(sReduced_fn),
                finfo = fitsinfo(sReduced_fn);
                S.ReducedKeys = finfo.PrimaryData.Keywords;
            end
            
            S.runnum = runnum;
            S.iter = iter;
            S.Reduced_fn = sReduced_fn;
            S.Rundir_fn = sRundir_fn;        

            % number of probes, and images per wavelength if probed
            % from lyotserver_IFS.py
            S.Nppair = (FitsGetKeywordVal(S.ImKeys,'NUMIM')+FitsGetKeywordVal(S.ImKeys,'PROBEF'));
            S.NumImProbe = 2*S.Nppair+1;
            % should check S.ImKeys, 'NOPAIRx' to see if wavelength x is probed

            % Norm Intensity from first, non-probed image, for each wave
            S.NofW = FitsGetKeywordVal(S.ImKeys,'NCOLOR');
            S.Nlamcorr = S.NofW; % not strictly corect, should come from reduced or config length(S.ilamcorr);
            S.imgindex = 1:S.NumImProbe:S.NofW*S.NumImProbe;
            
            % get actual wavelengths from the camera fits files
            for iwv = 1:S.NofW,
               fn = FitsGetKeywordVal(S.ImKeys, ['C' num2str(iwv-1) 'P0J0']) ;
               if ~exist(PathTranslator(fn),'file'),
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

            % default values from the class
            RminSc = S.RminSc;
            RmaxSc = S.RmaxSc;
            
            % parse options
            iRminSc = find(strcmpi('RminSc',varargin));
            if ~isempty(iRminSc),
                RminSc = varargin{iRminSc(end)+1};
            end
            iRmaxSc = find(strcmpi('RmaxSc',varargin));
            if ~isempty(iRmaxSc),
                RmaxSc = varargin{iRmaxSc(end)+1};
            end
                
            
            [x, y, X, Y, R] = CreateGrid(S.ImCubeUnProb{1}, 1./S.ppl0);
            bMaskSc = R >= RminSc & R <= RmaxSc;

            %if S.runnum == 603,
            %Sth = load('Throughput_20171003T122253.mat');
            cpath = regexp(S.Results_pn,'/','split');
            irr = find(strcmp('results',cpath));
            Sth = load(PathTranslator([strjoin({cpath{1:irr}},'/') '/Throughput_20171003T122253.mat']));

            dxyfudge = 0.9;
            
            rthpt = mean([Sth.src_dy(:).*dxyfudge./Sth.mmperlamD Sth.src_dx(:).*dxyfudge./Sth.mmperlamD],2);
            thpt = mean([Sth.thpt_dy(:) Sth.thpt_dx(:)],2);
            % normalize to clear region: lam/D > 8.5 & lam/D < 16.2
            thpt_norm = mean(thpt(abs(rthpt) >= 8.5 & abs(rthpt) < 16.2));
            thpt = thpt./ thpt_norm;
            % rthpt, thpt is now a lookup table to throughput for
            % cross-section
            % make a matrix
            Thpt = zeros(size(bMaskSc));
            Thpt(bMaskSc) = interp1(rthpt, thpt, R(bMaskSc));

%             for iwv = 1:S.NofW,
%                 ImCon{iwv} = S.ImCubeUnProb{iwv};
%                 iuse = R >= IWA & R <= OWA;
%                 ImCon{iwv}(~iuse) = 0.0;
%                 meanNI(iwv) = mean(ImCon{iwv}(iuse));
%                 ImCon{iwv}(iuse) = ImCon{iwv}(iuse) ./ interp1(rthpt, thpt, R(iuse));
%                 meancontrast(iwv) = mean(ImCon{iwv}(iuse));
%             end
            
            % note: don't like nonzeros(...) because pixels within the dark
            % hole that happen to have zero measured intesnity do not
            % contribute to the mean.
            Contrast.total_lam = zeros(S.NofW,1);
            for iwv = 1:S.NofW,
                Contrast.total_lam(iwv) = mean(nonzeros( S.ImCube(:,:,S.imgindex(iwv)).*S.bMask ));
                Contrast.score_lam(iwv) = mean(S.ImCubeUnProb{iwv}(bMaskSc));
                Contrast.contr_lam(iwv) = mean(S.ImCubeUnProb{iwv}(bMaskSc)./Thpt(bMaskSc));
                
            end
            Contrast.mean = mean(Contrast.contr_lam);

            % mean coherent and incoherent contrast
            
            if isempty(S.CohInt),
                S.ReadReducedCube;
            end

            % control region
            Contrast.inco_lam = zeros(S.NofW,1);
            Contrast.co_lam = zeros(S.NofW,1);
            for iwv = 1:S.Nlamcorr,
                bPampz = squeeze(S.bPampzero(iwv,:,:));
                %bMaskUse = ~bPampz & bMaskSc;
                bMaskUse = bMaskSc;
                %incinttmp = S.IncInt{iwv};
                %incinttmp(bPampz) = 0.0;
                %cohinttmp = S.CohInt{iwv};
                %cohinttmp(bPampz) = 0.0;
                %Contrast.inco_lam(iwv) = mean(nonzeros(incinttmp.*S.mdMask));
                %Contrast.co_lam(iwv)   = mean(nonzeros(cohinttmp.*S.mdMask));
                Contrast.inco_lam(iwv) = mean(S.IncInt{iwv}(bMaskUse)./Thpt(bMaskUse));
                Contrast.co_lam(iwv)   = mean(S.CohInt{iwv}(bMaskUse)./Thpt(bMaskUse));
            end
            Contrast.inco_mean = mean(Contrast.inco_lam);
            Contrast.co_mean = mean(Contrast.co_lam);
            
        end % GetContrast

        function [Excel, sheet] = ContrastReportExcel(S)
            
            if ~isequal(computer, 'PCWIN64'),
                Excel = [];
                sheet = [];
                warning('Must Use PCWIN To Open Excel');
                return
            end

            [Excel, Workbook, Sheets, sheet] = Exlopen;

            Exlsetval(sheet, {'A2','A6'}, [542 554 565 576 588].')
            
            RminSc = 6.3;
            RmaxSc = 7.3;
            col = 'B';
            Ctmp = S.GetContrast('RminSc',RminSc,'RmaxSc',RmaxSc);
            Exlsetval(sheet, {[col '2'],[col '6']}, Ctmp.contr_lam(:));
            Exlsetval(sheet, [col '7'], Ctmp.mean);
            
            RminSc = 7.3;
            RmaxSc = 18.5;
            col = 'C';
            Ctmp = S.GetContrast('RminSc',RminSc,'RmaxSc',RmaxSc);
            Exlsetval(sheet, {[col '2'],[col '6']}, Ctmp.contr_lam(:));
            Exlsetval(sheet, [col '7'], Ctmp.mean);
            
            RminSc = 18.5;
            RmaxSc = 19.5;
            col = 'D';
            Ctmp = S.GetContrast('RminSc',RminSc,'RmaxSc',RmaxSc);
            Exlsetval(sheet, {[col '2'],[col '6']}, Ctmp.contr_lam(:));
            Exlsetval(sheet, [col '7'], Ctmp.mean);
            
            RminSc = 6.3;
            RmaxSc = 19.5;
            col = 'E';
            Ctmp = S.GetContrast('RminSc',RminSc,'RmaxSc',RmaxSc);
            Exlsetval(sheet, {[col '2'],[col '6']}, Ctmp.contr_lam(:));
            Exlsetval(sheet, [col '7'], Ctmp.mean);

            
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
                isl = S.Nppair + 2*(1 + 1 + 1 + S.Nbscan + 1+1)-1;
                S.bPampzero(iwl,:,:) = logical(RedData(:,:,isl * S.Nlamcorr+iwl)); % + 1i*RedData(:,:,(isl + 1)*S.Nlamcorr+iwl);

                % mean coherent and incoherent contrast
                S.Contrast_inco(iwl) = mean(nonzeros(S.IncInt{iwl}.*S.mdMask));
                S.Contrast_co(iwl)   = mean(nonzeros(S.CohInt{iwl}.*S.mdMask));
                

                
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
                    ImPrPlus{iwv,ipr} = squeeze(S.ImCube(:,:,imnum));
                    ImPrMinus{iwv,ipr} = squeeze(S.ImCube(:,:,imnum+1));
                end
            end
            
            % Useful Derived Parameters                       
            
        end % ReadImageCube
        
        function S = ReadDMvCube(S)
            
            for dmnum = 1:2,
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
        
        function [hfig, ha] = DisplayImCubeUnProb(S, sOptions)
            
            if isempty(S.ImCube),
                S.ReadImageCube;
            end
            
            % default options and set requested options
            sOpt = struct(...
                'bPlotLog', false ...
                ,'drawRadii', [] ...
                ,'clim', [] ...
                ,'ilam', [] ...
                );
            if exist('sOptions','var') && isstruct(sOptions),
                fnames = fieldnames(sOptions);
                for ii = 1:length(fnames),
                    sOpt.(fnames{ii}) = sOptions.(fnames{ii});
                end
            end
            
            [x, y, X, Y, R] = CreateGrid(S.ImCubeUnProb{1}, 1./S.ppl0);
            xlim = 25*[-1 1]; ylim = xlim;
            
            if isempty(sOpt.ilam),
                sOpt.ilam = 1:S.NofW;
            end
            Nlam = length(sOpt.ilam);
            
            hfig = figure_mxn(1,Nlam);
            for iwv = 1:Nlam,
                iwvpl = sOpt.ilam(iwv);
                ha(iwv) = subplot(1,Nlam,iwv);
                if sOpt.bPlotLog,
                    imageschcit(x, y, log10(abs(S.ImCubeUnProb{iwvpl}))), axis image,
                    colorbartitle('log_{10}(Intensity)')
                else
                    imageschcit(x, y, S.ImCubeUnProb{iwvpl}), axis image,
                    colorbar
                end
               clim(iwv,:) = get(gca,'clim');
               
               set(gca,'xlim',xlim,'ylim',ylim)
               xlabel('\lambda/D'), ylabel('\lambda/D')
               title(['Iter #' num2str(S.iter) ', Wave ' num2str(S.NKTcenter(iwvpl)/S.NM) 'nm'])
                
               % overlay circles if requested
               if ~isempty(sOpt.drawRadii),
                   hold on
                   for irad = 1:length(sOpt.drawRadii),
                       draw_circle([0 0], 2*sOpt.drawRadii(irad), 1, 'r');
                   end
                   hold off
               end
               
            end % for each wavelength and subplot
            
            % set each image plot to the same clim
            % auto-clim, unless specific clim requested
            if isempty(sOpt.clim),
                cmaxsort = sort(clim(:),'descend');
                set(ha,'clim', [min(clim(:)) cmaxsort(2)]);
            else
                set(ha,'clim',sOpt.clim);
            end

            % radial plot
            Nr = 128;
            re = linspace(0, 25, Nr+1)';
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
            if ~isempty(sOpt.drawRadii),
                ylim = get(gca,'ylim');
                hold on
                for irad = 1:length(sOpt.drawRadii),
                    plot(sOpt.drawRadii(irad)*[1 1], ylim, '--r')
                end
                hold off
            end

            xlabel('Radius (\lambda/D)')
            ylabel('Average Normalized Intensity')
            legend(legstr{:},'location','north')
            title(['Iter #' num2str(S.iter)])
            
        end % DisplayImCubeUnProb(S)
        
        function [hfig, hax] = DisplayIncInt(S, varargin)
            % [hfig, hax] = DisplayIncInt(S, varargin)
            % options:
            %    'drawradii'
            %    'bLog'
            %    'clim'
            
            if isempty(S.IncInt),
                S.ReadReducedCube;
            end

            % options:
            iblog = strcmpi('blog',varargin);
            if any(iblog),
                ic = find(iblog);
                bLog = varargin{ic+1};
            else
                bLog = false;
            end
            
            % check option draw radii
            ibDrawRadii = strcmpi('drawradii',varargin);
            bDrawRadii = any(ibDrawRadii);
            if bDrawRadii,
                ic = find(ibDrawRadii);
                radii = varargin{ic+1};
            end
            
            iClim = strcmpi('clim',varargin);
            bClim = any(iClim);
            if bClim,
                climopt = varargin{find(iClim)+1};
            end
            
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
            
            hfig = figure_mxn(1,S.Nlamcorr);
            
            xlim = 25*[-1 1]; ylim = xlim;

            if bClim,
                clim = pFun(climopt);
            else
                clim = pClim([S.IncInt{:}]);
            end
            
            [x, y] = CreateGrid(S.IncInt{1}, 1./S.ppl0);
            % auto-scale
            %Agg = [
            for iwv = 1:S.Nlamcorr,
            
                hax(iwv) = subplot(1, S.Nlamcorr, iwv);
                him(iwv) = imageschcit(x, y, pFun(S.IncInt{iwv})); axis image
                xlabel('\lambda / D')
                ylabel('\lambda / D')
                %caxis(clim);
                colormap(gray)
                colorbartitle(cbartitle)
                title(['Inc Int, \lambda = ' num2str(S.NKTcenter(iwv)/S.NM) 'nm, it#' num2str(S.iter)])
            
            end
            set(hax,'xlim',xlim,'ylim',ylim,'clim',clim);
            
            if bDrawRadii,
                for iax = 1:length(hax),
                    axes(hax(iax));
                    hold on
                    for irad = 1:length(radii),
                        draw_circle([0 0], 2*radii(irad), 1, 'r');
                    end
                    hold off
                end % for each axes               
            end % if option drawradii
            
            % radial plot
            for iwv = 1:S.Nlamcorr,
                [fovrplot, IncIntrad] = RadialMean(x, y, S.IncInt{iwv}, 128);
                % remove negative inc int values so we can use log plot
                IncIntrad(IncIntrad < 1e-11) = 1e-11;
                qplot{1,iwv} = fovrplot;
                qplot{2,iwv} = IncIntrad;
                legstr{iwv} = [num2str(S.NKTcenter(iwv)/S.NM) 'nm'];
            end
            figure, semilogy(qplot{:}); 
            hax(end+1) = gca;
            grid on
            set(gca,'xlim',[0 1].*xlim)
            if bClim,
                set(gca,'ylim',clim)
            end
            xlabel('Radius (\lambda/D)')
            ylabel('Unmodulated (Norm. Int.)')
            title(['Inc Int, it#' num2str(S.iter)])
            
            if bDrawRadii,
                hold on
                for irad = 1:length(radii),
                    plot(radii(irad)*[1 1], get(gca,'ylim'), '--r')
                end
                hold off
            end
            
            % legend must go after everything is plotted
            legend(legstr{:},'Location','North')

        end % DisplayIncInt
                
        function hfig = DisplayProbeAmp(S)
            if isempty(S.ProbeAmp),
                S.ReadReducedCube;
            end
            
            xlimlamD = 22*[-1 1];
            [x, y] = CreateGrid(S.ProbeAmp{1,1}, 1./S.ppl0);
            
            % ProbeAmp is from reduced data cube, Primary HDU
            % arrange subplots
            Nwide = max([S.Nppair, S.Nlamcorr]);
            Nhigh = min([S.Nppair, S.Nlamcorr]);
            hfig = figure_mxn(Nhigh, Nwide);
            for ip = 1:S.Nppair,
                for iw = 1:S.Nlamcorr,
                    if Nwide == S.Nlamcorr,
                        ipl = (ip-1)*S.Nlamcorr+iw;
                    else
                        ipl = (iw-1)*S.Nppair+ip;
                    end
                    subplot(Nhigh, Nwide, ipl),
                    imageschcit(x, y, S.ProbeAmp{iw,ip}), axis image, colorbar,
                    set(gca,'xlim',xlimlamD,'ylim',xlimlamD)
                    title(['Iteration # ' num2str(S.iter) ' Wave #' num2str(iw) ' Probe # ' num2str(ip) ' Amplitude'])
                    
                end
            end

            % again, but only show points where bPampzero = false
            % i.e. pixels used by EFC
            hfig = figure_mxn(Nhigh, Nwide);
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
                    imageschcit(x, y, probeamptmp), axis image, colorbar,
                    set(gca,'xlim',xlimlamD,'ylim',xlimlamD)
                    title(['Iteration # ' num2str(S.iter) ' Wave #' num2str(iw) ' Probe # ' num2str(ip) ' Amplitude'])
                    
                end
            end

        end % DisplayProbeAmp
        
        function [pl, hfig] = DisplayOneProbeAmp(S, ip, drawradii)
            if isempty(S.ProbeAmp),
                S.ReadReducedCube;
            end
            if nargin < 3, drawradii = []; end
                
            % display one probe at all wavelengths and plot cross sections
            [x, y] = CreateGrid(S.ProbeAmp{1,1}, 1./S.ppl0);
            xlimlamD = 22*[-1 1];
            hfig = figure_mxn(1,S.Nlamcorr);
            for iw = 1:S.Nlamcorr,
                subplot(1,S.Nlamcorr,iw);
                imageschcit(x, y, (S.ProbeAmp{iw, ip}).^2), axis image, colorbartitle('Norm Intensity')
                set(gca,'xlim',xlimlamD,'ylim',xlimlamD)
                title(['Iteration # ' num2str(S.iter) ' Wave #' num2str(iw) ' Probe # ' num2str(ip) ' Amplitude'])
                
                if ~isempty(drawradii),
                    hold on
                    for ir = 1:length(drawradii),
                        draw_circle([0 0],2*drawradii(ir),1, 'r');
                    end
                    hold off
                end
                
                pl{1,iw} = [x(:) y(:)];
                pl{2,iw} = [S.ProbeAmp{iw, ip}(y==0, :).' S.ProbeAmp{iw, ip}(:, x==0)].^2;
                
            end 
            
            figure, semilogy(pl{:}), grid
            set(gca,'xlim',xlimlamD)
            xlabel('X, Y (\lambda/D)')
            ylabel('Norm Intensity')
            title(['Iter # ' num2str(S.iter) ', Probe # ' num2str(ip) ])
            
            
        end % DisplayOneProbeAmp
        
        function [hfig, ha] = DisplayProbeCube(S, iwvplot)
            % hfig = DisplayProbeCube(S, iwvplot)
            % ProbeCube is 2nd HDU of reduced data cube
            if isempty(S.ProbeModel),
                S.ReadProbeCube;
            end
            if ~exist('iwvplot','var') || isempty(iwvplot),
                warning('wavelength # not specified, plotting wave #1');
                iwvplot = 1;
            end
            
            [x, y] = CreateGrid(S.ProbeModel{1,1}, 1./S.ppl0);
            
            ip = 1;
            figure, plot(y, [abs(S.ProbeModel{iwvplot,ip}(:,x==0)), ...
                    S.ProbeMeasAmp{iwvplot,ip}(:,x==0)].^2), grid
            xlabel('Y (\lambda / D)')
            legend(['it ' num2str(S.iter) ' model'],['it ' num2str(S.iter) ' Measured'])
            title(['Cross Section, it ' num2str(S.iter) ', Wave #' num2str(iwvplot) ', Measure Probe #' num2str(ip)])
            
            hfig = figure_mxn(2,S.Nppair);
            for ip = 1:S.Nppair,
                ha(ip) = subplot(2,S.Nppair,ip);
                imageschcit(x,y,abs(S.ProbeModel{iwvplot,ip}).^2), axis image, colorbar
                set(gca,'xlim',25*[-1 1],'ylim',25*[-1 1])
                xlabel('\lambda / D'), ylabel('\lambda / D')
                title(['wave #' num2str(iwvplot) ', it ' num2str(S.iter) ', Model Probe #' num2str(ip)])
                
                ha(S.Nppair+ip) = subplot(2,S.Nppair,S.Nppair+ip);
                imageschcit(x,y,S.ProbeMeasAmp{iwvplot,ip}.^2), axis image, colorbar
                set(gca,'xlim',25*[-1 1],'ylim',25*[-1 1])
                xlabel('\lambda / D'), ylabel('\lambda / D')
                title(['wave #' num2str(iwvplot) ', it ' num2str(S.iter) ', Measure Probe #' num2str(ip)])
            end

        end % DisplayProbeCube
        
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
            Nfft = Nlam2 ./(S.NKTcenter(2)./S.NKTcenter);
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
            xlim = 25*[-1 1]; ylim = xlim;

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
            
            if isempty(S.DMvCube)
                S.ReadDMvCube;
            end
            
            switch length(varargin),
                case 0,
                    refDM1v = 0;
                    strRefDM1 = '';

                    refDM2v = 0;
                    strRefDM2 = '';
                    
                case 1,
                    if isa(varargin{1},'CRunData'),
                        Sref = varargin{1};
                        if isempty(Sref.DMvCube)
                            Sref.ReadDMvCube;
                        end
                        
                        refDM1v = squeeze(Sref.DMvCube{1}(:,:,1));
                        strRefDM1 = ['it#' num2str(Sref.iter)];
                        refDM2v = squeeze(Sref.DMvCube{2}(:,:,1));
                        strRefDM2 = ['it#' num2str(Sref.iter)];
                        
                    else, 
                        error('usage');
                    end
                    
                case 2,
                    [refDM1v_fn, refDM2v_fn] = deal(varargin{1:2});
                    
                    refDM1v = fitsread(PathTranslator(refDM1v_fn));
                    aatmp = regexp(refDM1v_fn, '/', 'split');
                    strRefDM1 = pwd2titlestr(aatmp{end});

                    
                                
                    refDM2v = fitsread(PathTranslator(refDM2v_fn));
                    aatmp = regexp(refDM2v_fn, '/', 'split');
                    strRefDM2 = pwd2titlestr(aatmp{end});


                otherwise %else
                    error('usage');
                    
            end % switch nargin

            hfig = figure_mxn(1,2);
            
            dDMv1 = squeeze(S.DMvCube{1}(:,:,1)) - refDM1v;
            rmsdDMv1 = rms(dDMv1(abs(dDMv1)>0));
            hax(1) = subplot(1,2,1);
            imageschcit(dDMv1)
            colorbartitle('Vmu')
            title(['Iter # ' num2str(S.iter) '; DM1; \DeltaVmu, ref = ' strRefDM1])
            
            dDMv2 = squeeze(S.DMvCube{2}(:,:,1)) - refDM2v;
            rmsdDMv2 = rms(dDMv2(abs(dDMv2)>0));
            hax(2) = subplot(1,2,2);
            imageschcit(dDMv2)
            colorbartitle('Vmu')
            title(['Iter # ' num2str(S.iter) '; DM2; \DeltaVmu, ref = ' strRefDM2])
            
            fprintf('rms dDMv1 = %.3f Vmu\n',rmsdDMv1);
            fprintf('rms dDMv2 = %.3f Vmu\n',rmsdDMv2);
            
        end % DisplayDMv
        
        %
        %         function DisplayRadialContrast(S)
        %
        %         end % DisplayRadialContrast
        
    end % methods
    
end % classdef