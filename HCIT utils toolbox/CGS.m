classdef CGS < handle
    % S = CGS(gsnum, bn)
    %
    % bn choices, or specify:
    %                     case 'dst'
    %                         bn = '/proj/dst/data/dB_PR/gsdst_';
    %                     case 'spc_disc'
    %                         bn = '/home/dmarx/HCIT/SPC_disc/gsspc_20171204/reduced/gsspc_';
    %                     case 'mcb_spc'
    %                         bn = '/home/dmarx/HCIT/MCB_SPC/phaseretrieval/reduced/gsomc_no';
    %                     case 'mcb_hlc'
    %                         bn = '/proj/mcb/data/dB_PR_Kern/gsomc_no00';
    %                     case 'ttb_hlc'
    %                         bn = '/proj/mcb/data/dB_PR_Kern/gsomc_no00';
    %
    % read reduced results from gs phase retrieval
    %
    % methods:
    %    [hfig, hax] = DisplayGS(S)
    %
    %    [hfig, hax] = DisplayGSrefGS(S, Sref)
    %          hfig = CheckOption('hfig', figure_mxn(2,2), varargin{:});
    %          usebMask = CheckOption('usebMask', true, varargin{:});
    %          phplot = CheckOption('phplot', 'angleE', varargin{:}); % or S.(phplot)
    %          dphclim = CheckOption('dphclim', [], varargin{:});
    %
    %    [Z, rz, pharesidual] = ZernikeFit(S, nz)
    %
    %    Scopy = Copy(S)
    %
    %    cc = AmpCorrMetric(S)
    %
    %    [hfig, hax] = DisplayAmpPlane(S, ipl)
    %
    % hax = DisplayAllPlanes(S, options)
    %      options:
    %         'image': 'meas' (default) or 'calc'
    %         'value': 'amp' (default) or 'intensity'
    %
    % properties (some of them):
    %    cAmpPlanes{1:Num Images}
    %       % (:,:,1) = measured amplitude
    %       % (:,:,2) = calculated amplitude
    %       % (:,:,3) = calculated phase

    properties
        
        gsnum 
        listPupImDir
        listSrcImDir
        bn 
        amp
        ph         % unwrapped phase
        phunwrap_bMask % mask used in unwrap (WFSC version)
        phw        % wrapped phase
        cAmpPlanes % input and est amp images each plane (cmp.fits)
        zAmpPlanes % camz for each amp image
        amp_keys
        bMask
        ampthresh
        phunwrap
        phw_ptt
        E
        x
        y
        X
        Y
        R
        T
        wavelength

        % Remap is for PIAA (upstream of PIAA tube is remapped to downstream)
        RemapRadialR2
        RemapRadialR1
        RemapRadialRmin % = 1.51*MM;
        RemapRadialRmax %= 14.2*MM;
        RemapRadialRpix  % = 145; % pixels pupil radius taken manually from bMask(:,x==0)
        Eremap
        phunwrapremap
        bMaskRemap
        
    end % properties
    
    methods
        
        function S = CGS(gsnum, bn, varargin)
            % S = CGS(gsnum, bn)
            %
            % example:
            %
            
            if nargin == 0,
                return
            end
            
            U = CConstants;
            
            if ~exist('bn','var') || isempty(bn),
                %bn = '/home/dmarx/HCIT/DST/phaseretrieval_20180605/reduced/gsdst_';
                bn = '/proj/dst/data/dB_PR/gsdst_';
            else
                switch lower(bn),
                    case 'dst'
                        bn = '/proj/dst/data/dB_PR/gsdst_';
                        if gsnum <= 284,
                            year = '2019';
                            uname = 'bseo';
                        elseif gsnum <= 334,
                            year = '2019';
                            uname = 'bseo';
                        else
                            year = '2019';
                            uname = 'dmarx';
                        end
                        
                        S.listPupImDir = dir(PathTranslator(...
                            ['/proj/piaa-data/Data/' year '-*-*/' uname '/gsdst_p_' num2str(gsnum,'%04d') '/*.fits'] ...
                            ));
                        S.listSrcImDir = dir(PathTranslator(...
                            ['/proj/piaa-data/Data/' year '-*-*/' uname '/gsdst_s_' num2str(gsnum,'%04d') '/*.fits'] ...
                            ));
                       
                    case 'spc_disc'
                        bn = '/home/dmarx/HCIT/SPC_disc/gsspc_20171204/reduced/gsspc_';
                    case 'mcb_spc'
                        bn = '/home/dmarx/HCIT/MCB_SPC/phaseretrieval/reduced/gsomc_no';
                        % get dir listing of raw camera images
                        S.listPupImDir = dir(PathTranslator(...
                            ['/proj/piaa-data/Data/2019-*-*/dmarx/gsomc_p_' num2str(gsnum) '/*.fits']...
                            ));
                        S.listSrcImDir = dir(PathTranslator(...
                            ['/proj/piaa-data/Data/2019-*-*/dmarx/gsomc_s_' num2str(gsnum) '/*.fits']...
                            ));
                    case 'mcb_hlc'
                        bn = ['/proj/mcb/data/dB_PR_Kern/gsomc_no' num2str(gsnum,'%05d')];

                        % get dir listing of raw camera images
                        if gsnum >= 843,
                            year = '2019';
                            gsbn = 'gsttb';
                        elseif gsnum >= 806,
                            year = '2018';
                            gsbn = 'gsomc';
                        else
                            error('which year is this?');
                        end
                        S.listPupImDir = dir(PathTranslator(...
                            ['/proj/piaa-data/Data/' year '-*-*/bseo/' gsbn '_p_' num2str(gsnum,'%05d') '/*.fits']...
                            ));
                        S.listSrcImDir = dir(PathTranslator(...
                            ['/proj/piaa-data/Data/' year '-*-*/bseo/' gsbn '_s_' num2str(gsnum,'%05d') '/*.fits']...
                            ));

                        % for EMCCD, starting, 2020-11-04, gsnum 1369
                        wavelength_kwd = 'lam';

                    case 'mcb_alllens'
                        % optional trial name because we might be testing different ways
                        % to process one gsnum
                        trialname = '';
                        if ~isempty(varargin),  trialname = [varargin{1} '_']; end
                        bn = ['/home/dmarx/WFIRST/PR_lead/all_lens_PR/all_lens_20201203/reduced/prout_' trialname num2str(gsnum)];
                                                
                        S.listPupImDir = dir(PathTranslator(...
                            ['/proj/mcb/data/excam/2021-*-*/gsomc_p_' num2str(gsnum,'%d') '/*.fits']...
                            ));
                        S.listSrcImDir = dir(PathTranslator(...
                            ['/proj/mcb/data/excam/2021-*-*/gsomc_s_' num2str(gsnum,'%d') '/*.fits']...
                            ));

                        wavelength_kwd = 'lam';
                        
                    case 'mcb_twolens'
                        % optional trial name because we might be testing different ways
                        % to process one gsnum
                        trialname = '';
                        if ~isempty(varargin),  trialname = [varargin{1} '_']; end
                        bn = ['/home/dmarx/WFIRST/PR_lead/all_lens_PR/two_lens_20201109/reduced/prout_' trialname num2str(gsnum)];
                                                
                        S.listPupImDir = dir(PathTranslator(...
                            ['/proj/mcb/data/excam/2021-*-*/gsomc_p_' num2str(gsnum,'%d') '/*.fits']...
                            ));
                        S.listSrcImDir = dir(PathTranslator(...
                            ['/proj/mcb/data/excam/2021-*-*/gsomc_s_' num2str(gsnum,'%d') '/*.fits']...
                            ));

                        wavelength_kwd = 'lam';

                    case 'ttb_hlc'
                        bn = ['/proj/mcb/data/dB_PR_Kern/gsomc_no' num2str(gsnum,'%05d')];
                        % get dir listing of raw camera images
                        if gsnum >= 843,
                            year = '2019';
                        else
                            error('which year is this?');
                        end
                        S.listPupImDir = dir(PathTranslator(...
                            ['/proj/piaa-data/Data/' year '-*-*/bseo/gsttb_p_' num2str(gsnum,'%05d') '/*.fits']...
                            ));
                        S.listSrcImDir = dir(PathTranslator(...
                            ['/proj/piaa-data/Data/' year '-*-*/bseo/gsttb_s_' num2str(gsnum,'%05d') '/*.fits']...
                            ));

                    case 'piaacmc'
                        % optional trial name because we might be testing different ways
                        % to process one gsnum
                        trialname = '';
                        if ~isempty(varargin),  trialname = [varargin{1} '_']; end

                        bn = ['/proj/piaacmc/phaseretrieval/reduced/piaa_' trialname num2str(gsnum,'%03d')];
                        
                        % get dir listing of raw camera images                        
                        S.listPupImDir = dir(PathTranslator(...
                            ['/proj/piaacmc/scicam/*/gspiaa_p_' num2str(gsnum,'%04d') '/piaa*.fits']...
                            ));
                        S.listSrcImDir = dir(PathTranslator(...
                            ['/proj/piaacmc/scicam/*/gspiaa_s_' num2str(gsnum,'%04d') '/piaa*.fits']...
                            ));
                        
                    otherwise
                        % do nothing, let bn = bn
                        bn = [bn num2str(gsnum,'%d')];
                        %disp(bn);
                end % switch bn
            end % if ~exist('bn','var') || isempty(bn),
            
            %disp(['opening: ' PathTranslator([bn num2str(gsnum,'%03d') 'amp.fits'])]);
            try
                %ampinfo = fitsinfo(PathTranslator([bn num2str(gsnum,'%03d') 'amp.fits']));
                ampinfo = fitsinfo(PathTranslator([bn 'amp.fits']));
            catch
                fprintf('failed to open fits file: \n%s\n',PathTranslator([bn 'amp.fits']));
            end
            
 
            S.gsnum = gsnum;
            S.bn = bn;
            S.amp = fitsread(PathTranslator([bn 'amp.fits']));
            S.ph = fitsread(PathTranslator([bn 'ph.fits']));  % unwrapped phase
            S.phw = fitsread(PathTranslator([bn 'phwrap.fits'])); % = angle(eref), wrapped phase
            S.amp_keys = ampinfo.PrimaryData.Keywords;
            %S.wavelength = FitsGetKeywordVal(S.amp_keys, wavelength_kwd);
                       
            % phase unwrap in WFSC puts the mask used for PR in the second
            % hdu (starting Feb 2021)
            finfo = fitsinfo(PathTranslator([bn 'ph.fits']));
            if length(finfo.Contents) > 1,
                fn = [bn 'ph.fits'];
                S.phunwrap_bMask = logical(fitsread(PathTranslator(fn),'image'));
            end
            
            
            % S.bMask, S.ampthresh
            [sResult, S.bMask] = AutoMetric(S.amp, [], struct('AutoThreshold_Nbins',32));
            S.ampthresh = sResult.thresh;
            % morphological opening to eliminate outlier 'salt' noise
            try
                S.bMask = imclose(S.bMask, strel('disk',3));
                S.bMask = imopen(S.bMask, strel('disk',3));
            catch
                warning('maybe no image processing license, skipping imclose, imopen for bMask');
            end

            % check that mask is reasonable
            [B,L,N,A] = bwboundaries(S.bMask, 'noholes');
            if N > 10,
                warning(['pupil mask has ' num2str(N) ' regions']);
                %keyboard;
            end

            % 
            %S = AdjustUnwrapRegionPiston(S);
            
            % S.phw_ptt
            % use FFT to remove large amounts of PTT (integer pixels in FFT space
            % then use zernikes to remove remaining PTT
            [S.x, S.y, S.X, S.Y, S.R, S.T] = CreateGrid(S.amp);
            S.RemovePTTfft; % creates first estimate of S.phw_ptt
            S.phw_ptt = RemovePTTZ(S.phw_ptt, S.bMask);
            %             S.phunwrap = RemovePTTZ(S.phunwrap, S.bMask);
            %             S.phw_ptt  = mod2pi(S.phunwrap);

            %             % unwrap phase using better unwrap routine, but requires mask
            %             %phw = S.phw;
            %             phw = S.phw_ptt + 1 - 1;
            %             phw(~S.bMask) = NaN;
            %             S.phunwrap = unwrap_phase(phw);
            %             S.phunwrap(~S.bMask) = 0;
            S.phunwrap = S.phw_ptt;

            % S.E
            S.E = S.amp .* exp(1i*S.phw_ptt);
            %S.E = S.amp .* exp(1i*S.phw);
        
            % load the radial mapping
            % should be part of the bn switch
            % see email from Dan Sirbu "RE double check r2 r1 definition"
            try
                rm2_rm1 = load(PathTranslator('/proj/piaacmc/phaseretrieval/2019-10-16-nutekPiaaRemappingCoords/remapping.txt'));
                S.RemapRadialR2 = rm2_rm1(:,1);
                S.RemapRadialR1 = rm2_rm1(:,2);
                S.RemapRadialRmin = 1.51*U.MM;
                S.RemapRadialRmax = 14.2*U.MM;
                % 145 pixels pupil radius taken manually from bMask(:,x==0)
                % 15.0mm ray trace mag * measured 0.5* 46.3mm beam diameter at
                % pupil
                S.RemapRadialRpix = 192.5;
            catch
                warning('could not load remapping.txt');
            end
            
        end % CGS instantiator
        
        function Scopy = Copy(S)
            Scopy = CGS;
            
            listFnames = fieldnames(S);
            for ii = 1:length(listFnames),
                Scopy.(listFnames{ii}) = S.(listFnames{ii});
            end

        end
        
        function Crop(S, wpix)
                       
            % crop all the image arrays to something bigger than the mask
            if nargin == 0 || isempty(wpix),
                wpix = 1.5*max(S.R(S.bMask));
            end
            
            % guarantee that image size is even number
            wc = ceil(wpix/2);
            
            [S.amp, S.bMask] = CropImage(S.amp, S.bMask, [0 0], 2*wc, 2*wc);
            S.ph  = CropImage(S.ph, [], [0 0], 2*wc, 2*wc);
            S.phw = CropImage(S.phw, [], [0 0], 2*wc, 2*wc);
            S.phw_ptt = CropImage(S.phw_ptt, [], [0 0], 2*wc, 2*wc);
            S.phunwrap = CropImage(S.phunwrap, [], [0 0], 2*wc, 2*wc);            
            S.E   = CropImage(S.E, [], [0 0], 2*wc, 2*wc);            
            
            [S.x, S.y, S.X, S.Y, S.R, S.T] = CreateGrid(S.amp);

        end % Crop
        
        function ReadAmpImages(S)
            % read the cmp.fits file
            cmp_fn = [S.bn 'cmp.fits'];
            
            if ispc,

                [~, fntmp, ~] = fileparts(S.bn);
                cmp_fn_local = ['C:\Users\dmarx\HCITdatatemp\' fntmp 'cmp.fits'];
                copyfile(PathTranslator(cmp_fn), cmp_fn_local);
                cmp_fn = cmp_fn_local;

            end
            
            finfo = fitsinfo(PathTranslator(cmp_fn));            
            NIm = length(finfo.Contents);

            S.cAmpPlanes = cell(1,NIm);
            S.zAmpPlanes = zeros(NIm,1);
            
            % each hdu is N x N x 3
            % (:,:,1) = measured amplitude
            % (:,:,2) = calculated amplitude
            % (:,:,3) = calculated phase
            
            S.cAmpPlanes{1} = fitsread(PathTranslator(cmp_fn));
            S.zAmpPlanes(1) = FitsGetKeywordVal(finfo.PrimaryData.Keywords, 'Z');
            for ii = 2:NIm,
                S.cAmpPlanes{ii} = fitsread(PathTranslator(cmp_fn),'image',ii-1);
                S.zAmpPlanes(ii) = FitsGetKeywordVal(finfo.Image(ii-1).Keywords, 'Z');
            end
            
        end % ReadAmpImages
        
        function rmsP = rmsPha(S)
            
            rmsP = rms(S.phunwrap(S.bMask));
            
        end % rmsPha
        
        function rmsP = rmsPhaDiff(S, Sref)

            % intersection of mask
            bMaskTmp = S.bMask & Sref.bMask;
            
            % S.phunwrap(S.bMask) is already zero-mean
            % if S.bMask is very different from Sref.bMask, then
            % difference might not be zero mean            
            phdiff = S.phunwrap(bMaskTmp) - Sref.phunwrap(bMaskTmp);
            rmsP = rms(phdiff - mean(phdiff));
            
        end % rmsPha

        function [hfig, hax] = DisplayGS(S, varargin)
            % [hfig, hax] = DisplayGS(S, varargin)
            %
            % pMask = CheckOption('pMask', S.bMask, varargin{:});
            % xylim = CheckOption('xylim', 1.1*max(S.R(S.bMask)), varargin{:});
            % climph = CheckOption('climph', [], varargin{:});
            % phplot = CheckOption('phplot', 'angleE', varargin{:}); %
            %     other choices = 'phw_ptt', 'ph', 'phw', 'phunwrap'
            % ampplot = CheckOption('ampplot', 'absE', varargin{:});
            %     other choices 'amp'
            % stitle = CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});

            pMask = CheckOption('pMask', S.bMask, varargin{:});
            xylim = CheckOption('xylim', 1.1*max(S.R(S.bMask)), varargin{:});
            climph = CheckOption('climph', [], varargin{:});
            phplot = CheckOption('phplot', 'angleE', varargin{:}); % or S.(phplot)
            ampplot = CheckOption('ampplot', 'absE', varargin{:});
            stitle = CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});
            
            %hfig = figure;
            %hax = imagescampphase(S.E, x, y, ['gsnum ' num2str(S.gsnum)]);

            %ampplot = CheckOption('ampplot', 'absE', varargin{:});
            switch lower(ampplot)
                case 'abse'
                    amp = abs(S.E); % default, masked
                otherwise
                    amp = S.(ampplot); % S.amp; not masked, right from fits file
            end
                   
            hfig = figure_mxn(1,2);
            hax(1) = subplot(1,2,1);
            imageschcit(S.x, S.y, abs(S.E))
            colorbartitle('Amplitude')
            set(gca,'xlim',xylim*[-1 1],'ylim',xylim*[-1 1])
            title(stitle)

            % bMask is only for phase plot
            % check options for what to plot
            switch phplot
                case 'angleE'
                    ph = angle(S.E);
                otherwise
                    ph = S.(phplot);
            end

            
            if isempty(pMask), pMask = ones(size(S.E)); end
            hax(2) = subplot(1,2,2);
            imageschcit(S.x, S.y, pMask.*ph)
            colorbartitle('Phase (rad)')
            set(gca,'xlim',xylim*[-1 1],'ylim',xylim*[-1 1])
            if isempty(climph),
                climph = AutoClim(ph, 'symmetric', true);
            end
            set(gca,'clim',climph)
            title(['rms\phi = ' num2str(S.rmsPha,'%.3f') 'rad'])
            
            
            
        end % DisplayGS

        function [hfig, hax, dphaResult] = DisplayGSrefGS(S, Sref, varargin)
            % [hfig, hax] = DisplayGSrefGS(S, Sref, options)
            %
            % options:
            %    ('hfig', figure_mxn(2,2), varargin{:});
            %    ('usebMask', true, varargin{:});
            %    ('removeDefocus', false, varargin{:});
            %    ('doRegister', false, varargin{:}); % (false), true
            %    ('phplot', 'angleE', (default) 'phw_ptt', 'phunwrap', S.(phplot)
            %    ('xylim', 1.1*max(S.R(S.bMask)), varargin{:});
            %    ('dphclim', [], varargin{:});
            
            % parse options
            hfig = CheckOption('hfig', figure_mxn(2,2), varargin{:});
            usebMask = CheckOption('usebMask', true, varargin{:});
            removeDefocus = CheckOption('removeDefocus', false, varargin{:});
            doRegister = CheckOption('doRegister', false, varargin{:});
            xylim = CheckOption('xylim', [], varargin{:});
            phplot = CheckOption('phplot', 'angleE', varargin{:}); % S.(phplot)
            climdph = CheckOption('dphclim', [], varargin{:});
            climph = CheckOption('climph', [], varargin{:});
            
            switch phplot
                case 'angleE'
                    funPhPl = @(S) angle(S.E);
                otherwise
                    funPhPl = @(S) S.(phplot);
            end
            
            figure(hfig);

            % determine plot width
            if isempty(xylim),
                xylim = 1.1*max(S.R(S.bMask));
                xylim = 5*ceil(xylim/5.0); % to the nearest multiple of 5
            end
            
            hax(1) = subplot(2,2,1);
            imageschcit(S.x, S.y, abs(S.E))
            colorbartitle('Amplitude')
            set(gca,'xlim',xylim*[-1 1],'ylim',xylim*[-1 1])
            title(['gsnum ' num2str(S.gsnum)])
            
            hax(2) = subplot(2,2,2);
            imageschcit(S.x, S.y, funPhPl(S))
            colorbartitle('Phase (rad)')
            set(gca,'xlim',xylim*[-1 1],'ylim',xylim*[-1 1])
            title(['gsnum ' num2str(S.gsnum)])
            if ~isempty(climph), set(gca,'clim', climph), end
            
            % if size(Sref.amp) ~= size(S.amp), make Sref same size as S
            if ~isequal(size(Sref.amp), size(S.amp))
                if all(size(Sref.amp) < size(S.amp)),
                    % pad Sref
                    warning('size does not match, padding Sref to match S');
                    Sreftmp = struct(...
                        'phplot', PadImArray(funPhPl(Sref), size(funPhPl(S))) ...
                        ,'amp', PadImArray(Sref.amp, size(S.amp)) ...
                        );                    
                elseif all(size(Sref.amp) > size(S.amp)),
                    % crop Sref
                    warning('size does not match, cropping Sref to match S');
                    [nr, nc] = size(S.amp);
                    Sreftmp = struct(...
                        'phplot', CropImage(funPhPl(Sref), [], [0 0], nc, nr) ...
                        ,'amp', CropImage(Sref.amp, [nr nc]) ...
                        );
                else
                    % not square? something is wrong
                    error(['sizes do not match, size(Sref.amp) = ' num2str(size(Sref.amp))]);
                end
                
            else % same size
                Sreftmp = struct(...
                    'phplot', funPhPl(Sref) ...
                    ,'amp', Sref.amp ...
                    );
                
            end
            
            % match COM, translate Sreftmp to match S
            if doRegister,
                [xt, yt, Xt, Yt] = CreateGrid(S.bMask);
                comS = calcCOM(xt, yt, S.bMask);
                [~, bMaskref] = AutoMetric(Sreftmp.amp);
                comRef = calcCOM(xt, yt, bMaskref);
                xyshift = comRef - comS;
                
                % % fft linear shift to translate
                % fDoShift = @(A) fftshift(ifft2(ifftshift(ifftshift(fft2(fftshift(A))).*exp(-1j*pi*(xyshift(1)*Xt/xt(1) + xyshift(2)*Yt/yt(1))))));
                % Etmp = fDoShift(Sreftmp.amp .* exp(1j*Sreftmp.phplot));
                % Sreftmp.phplot = angle(Etmp);
                % Sreftmp.amp    = abs(Etmp);
                
                % integer shift
                fDoShift = @(A, dxy) circshift( circshift(A, -dxy(2)), -dxy(1), 2);
                Sreftmp.phplot = fDoShift(Sreftmp.phplot, round(xyshift));
                Sreftmp.amp    = fDoShift(Sreftmp.amp,    round(xyshift));
                
            end
            
            % plot amplitude difference
            hax(3) = subplot(2,2,3);
            imageschcit(S.x, S.y, S.amp - Sreftmp.amp)
            colorbartitle('\Delta Amplitude')
            set(gca,'xlim',xylim*[-1 1],'ylim',xylim*[-1 1])
            title(['gsnum ' num2str(S.gsnum) ' Amp -  Ref gsnum ' num2str(Sref.gsnum) ' Amp'])
            
            
            % optional remove defocus difference
            if removeDefocus
                %                 [ZZ, rz, dphatmp] = S.ZernikeFit(1:4);
                %                 [ZZ, rz, dpharef] = Sref.ZernikeFit(1:4);

                Ediff = S.E .* Sreftmp.amp.*exp(-1j*Sreftmp.phplot);
                [ZZ, phaimg, dpha, sOptions] = ZernikeAnalysis(Ediff, 'modes', 1:4, 'bMask', S.bMask);
                dphaResult = struct(...
                    'ZZ', ZZ ...
                    ,'phaimg', phaimg ...
                    ,'dpha', dpha ...
                    ,'sOptions', sOptions ...
                    );
            else
                dpha = mod2pi(funPhPl(S) - Sreftmp.phplot);
                dphaResult = struct('dpha', dpha);
            end
            
            hax(4) = subplot(2,2,4);            
            if usebMask,
                dpha = S.bMask .* dpha;
            end
            him = imageschcit(S.x, S.y, dpha);
            colorbartitle('Phase (rad)')
            set(gca,'xlim',xylim*[-1 1],'ylim',xylim*[-1 1])
            if isempty(climdph), set(gca,'clim',AutoClim(dpha,'symmetric',true,'pctscale',100))
            else set(gca,'clim',climdph)
            end
            title(['gsnum ' num2str(S.gsnum) ' Ref gsnum ' num2str(Sref.gsnum) ', rms \Delta = ' num2str(rms(dpha(S.bMask)),'%.3f') 'rad'])
            
            
            
        end % DisplayGSrefGS

        function phw_ptt = RemovePTTfft(S)
            Zs = exp(1i*S.phw);

            ZZ = fftshift(fft2(fftshift(S.bMask.*Zs)));
            
            %figure, imagescampphase(ZZ, x, y)
            
            [fm, xm, ym] = findpeak2(abs(ZZ));
            apha = angle(ZZ(round(ym),round(xm)));
            
            ZZs = circshift( circshift(ZZ, -S.y(round(ym)), 1), -S.x(round(xm)), 2);
            %figure, imageschcit(S.x, S.y, abs(ZZs))
            
            % remove piston
            Zss = exp(-1i*apha).*fftshift(ifft2(fftshift(ZZs)));
            %figure, imagescampphase(Sspin.amp .* Zspinss)

            phw_ptt = angle(Zss);
            S.phw_ptt = phw_ptt;
            
        end % RemovePTTfft
        
        function [ZZ, rz, pharesidual] = ZernikeFit(S, nz, varargin)
            % [ZZ, rz, pharesidual] = ZernikeFit(S, nz, varargin)
            %
            % zernike fit using bMask pixes
            %
            % bDisplay = CheckOption('display', true, varargin{:});
            % titlestr = CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});

            bDisplay = CheckOption('display', true, varargin{:});
            titlestr = CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});
            xylim = CheckOption('xylim', [], varargin{:});
            phresclim = CheckOption('phresclim', [], varargin{:});
            
            % fit should always include piston, tip, tilt, even if not
            % included in requested nz
            nzfit = nz(:);
            if ~any(nzfit == 3), nzfit = [3; nzfit]; end
            if ~any(nzfit == 2), nzfit = [2; nzfit]; end
            if ~any(nzfit == 1), nzfit = [1; nzfit]; end
            
            %
            rz = max(S.R(S.bMask));
            ZZ = zernikefit(S.X(S.bMask), S.Y(S.bMask), S.phw_ptt(S.bMask), nzfit, rz, 'noll');
            phwfit = nan(size(S.X));
            phwfit(S.bMask) = zernikeval(ZZ, S.X(S.bMask), S.Y(S.bMask), rz, 'noll', 'nz', nzfit);
            pharesidual = nan(size(S.X));
            pharesidual(S.bMask) = mod2pi(S.phw_ptt(S.bMask) - phwfit(S.bMask));

            
            %%%%%% plot results
            if bDisplay,

                % determine plot width
                if isempty(xylim),
                    xylim = 1.1*max(S.R(S.bMask));
                    xylim = 5*ceil(xylim/5.0); % to the nearest multiple of 5
                    xylim = xylim*[-1 1];
                end

                
                figure_mxn(2,3)

                subplot(2,3,1), imageschcit(S.x, S.y, abs(S.E))
                title([titlestr '; Amplitude'])
                set(gca,'xlim',xylim,'ylim',xylim);

                % phase
                subplot(2,3,2), imageschcit(S.x, S.y, mod2pi(S.phunwrap).*S.bMask)
                colorbartitle('Phase (rad)')
                set(gca,'clim',pi*[-1 1])
                title([titlestr '; Phase rms\phi = ' num2str(S.rmsPha,'%.3f')])               
                set(gca,'xlim',xylim,'ylim',xylim);
                
                % residual phase
                subplot(2,3,3), imageschcit(S.x, S.y, pharesidual),
                colorbartitle('Phase (rad)')
                rmse = rms(pharesidual(S.bMask));
                title(['Residual Fit, rms error = ' num2str(rmse,'%.3f') 'rad'])
                set(gca,'xlim',xylim,'ylim',xylim);
                if ~isempty(phresclim), set(gca,'clim',phresclim); end
                
                % fit phase
                subplot(2,3,4), imageschcit(S.x, S.y, phwfit), 
                colorbartitle('Phase (rad)'), title('Fit Phase')
                set(gca,'xlim',xylim,'ylim',xylim);                
                
                % bar graph zernike order, don't plot piston
                subplot(2,3,5:6)
                [Zplot, NZplot] = filterdata(nzfit > 3, ZZ, nzfit);
                hh = bar(NZplot, Zplot); grid
                %                 set(gca,'XTick', 1:length(nz(2:end)))
                %                 set(gca, 'XTickLabel', num2str(nz(2:end)))

                ylabel('Zernike Coeff (rms rad)')
                xlabel('Zernike # (Noll Order)')
                set(gca,'ylim',2*[-1 1])

                
            end
            
        end
        
        function [cc, MF] = AmpCorrMetric(S, varargin)
            % [cc, MF] = AmpCorrMetric(S)
            %
            % cc(:,1) = amplitude correlation
            % cc(:,2) = intensity correlation
            %
            % MF = length(cc(:,1)) - sum(cc(:,1)) % add more options?
            
            if isempty(S.cAmpPlanes),
                S.ReadAmpImages;
            end
            
            CCor = @(a,b) a(:)'*b(:)./sqrt( (a(:)'*a(:)) * (b(:)'*b(:)) );
            CCorq = @(q) CCor(q(:,:,1),q(:,:,2));
            CCint = @(q) CCor(abs(q(:,:,1)).^2, abs(q(:,:,2)).^2);
            
            NIm = length(S.cAmpPlanes);
            cc = zeros(NIm,2);

            for ii = 1:NIm,
                cc(ii,1) = CCorq(S.cAmpPlanes{ii});
                cc(ii,2) = CCint(S.cAmpPlanes{ii});
            end

            % this is the merit function used to parameter searching
            MF = length(cc(:,1)) - sum(cc(:,1));
            
        end % AmpCorrMetric
        
        function [hfig, hax] = DisplayAmpPlane(S, ipl, varargin)
            % [hfig, hax] = DisplayAmpPlane(S, list_ipl)
            %
            % compare camera and estimated
            % options:
            %   Option('value', 'amp' (default), 'intensity'
            %   Option('clim', [], varargin{:});
            %   Option('scale', 'log', varargin{:});
            %   Option('xylim', [], varargin{:});
            
            if isempty(S.cAmpPlanes), S.ReadAmpImages; end
            
            plAmpOrInt = CheckOption('value', 'amp', varargin{:});
            clim = CheckOption('clim', [], varargin{:});
            scale = CheckOption('scale', 'log', varargin{:});
            xylim = CheckOption('xylim', [], varargin{:});
            
            % each hdu is N x N x 3
            % (:,:,1) = measured amplitude
            % (:,:,2) = calculated amplitude
            % (:,:,3) = calculated phase
            
            % scale
            switch scale
                case 'linear'
                    fScale = @(a) a;
                case {'log', 'log10'},
                    fScale = @(a) real(log10(a));
                otherwise
                    error('unknown scale');
            end
            
            % choose amp or intensity
            if strcmp(plAmpOrInt, 'amp')
                funPlot = @(A, icc) fScale(squeeze(A(:,:,icc)));
            elseif strcmp(plAmpOrInt, 'intensity')
                funPlot = @(A, icc) fScale(squeeze(A(:,:,icc)).^2);
            end
            
            Nplanes = length(ipl);
            hfig = figure_mxn(2, Nplanes);
            hax = zeros(2, Nplanes);
            for iipl = 1:Nplanes,
                hax(1,iipl) = subplot(2,Nplanes,iipl);
                imageschcit(funPlot(S.cAmpPlanes{ipl(iipl)},1)), title(['Camera, ipl = ' num2str(ipl(iipl))])
                colorbar
                hax(2,iipl) = subplot(2,Nplanes,Nplanes + iipl);
                imageschcit(funPlot(S.cAmpPlanes{ipl(iipl)},2)), title(['Estimated, ipl = ' num2str(ipl(iipl))])
                colorbar
            end
            
            cclim = get(hax,'clim');
            if isempty(clim),
                set(hax,'clim', [max([cclim{1}(1) cclim{2}(1)]) max([cclim{:}])]);
            else
                set(hax,'clim',clim)
            end
            
            if ~isempty(xylim)
                set(hax, 'xlim', xylim, 'ylim', xylim)
            end
            
        end % DisplayAmpPlane
        
        function [hax, Imgs] = DisplayAllPlanes(S, varargin)
            % hax = DisplayAllPlanes(S, options)
            % options:
            %   'image': 'meas' (default) or 'calc'
            %   'value': 'amp' (default) or 'intensity' or 'phase' (calc only)
            %   'blog':  false (default) or true
            
            % each hdu is N x N x 3
            % (:,:,1) = measured amplitude
            % (:,:,2) = calculated amplitude
            % (:,:,3) = calculated phase

            U = CConstants;
            
            if isempty(S.cAmpPlanes), S.ReadAmpImages; end

            % parse options
            plMeasOrCalc = CheckOption('image', 'meas', varargin{:});
            plAmpOrInt = CheckOption('value', 'amp', varargin{:});
            bLog = CheckOption('blog', false, varargin{:});

            % choose measured or calculated
            if strcmp(plMeasOrCalc, 'meas'),
                funCamp = @(A) squeeze(A(:,:,1));
                strLabel = 'Measured';
            elseif strcmp(plMeasOrCalc, 'calc')
                funCamp = @(A) squeeze(A(:,:,2));
                strLabel = 'Calculated';
            else
                error(['Unknown Option image: ' plMeasOrCalc]);
            end
            strLabel = ['gsnum ' num2str(S.gsnum) ', ' strLabel];
            
            % choose amp or intensity
            if strcmp(plAmpOrInt, 'amp')
                funPlot = @(A) funCamp(A);
            elseif strcmp(plAmpOrInt, 'intensity')
                funPlot = @(A) funCamp(A).^2;
            elseif strcmp(plAmpOrInt, 'phase')
                funPlot = @(A) squeeze(A(:,:,3));
            else
                error(['Unknown Option value: ' plAmpOrInt]);
            end

            
            Ni = length(S.cAmpPlanes);
            Nc = ceil(Ni/2);
            Nr = 2;
            figure_mxn(Nr, Nc);
            hax = zeros(Ni,1);
            
            for ii = 1:Ni,
                hax(ii) = subplot(Nr, Nc, ii);
                
                Im = funPlot(S.cAmpPlanes{ii});

                if bLog,
                    Imlog = log10(Im);
                    imageschcit(Imlog)
                    set(gca,'clim',[-4 0] + max(Imlog(:)))
                    
                else
                    imageschcit(Im)
                end
                
                title(['Z = ' num2str(S.zAmpPlanes(ii)/U.MM,'%.1f')])

                Imgs{ii} = Im;
                
            end % for 

            % add the label
            w = 0.2;
            ha = annotation('textbox', [0.5-0.5*w 0.5 w 0.01] ...
                ,'String', strLabel ...
                ,'HorizontalAlignment','center' ...
                ,'VerticalAlignment','middle' ...
                ,'FitBoxToText','on','LineStyle','-' ...
                ,'EdgeColor','r','LineWidth',2 ...
                ,'FontSize', 24, 'FontWeight', 'bold' ...
                ,'Color','k' ...
                );
            
        end % DisplayAllPlanes
        
        function [r, Ir, hax] = DisplayRadialIntensity(S, varargin)
            % display mean intensity vs radius
            
            [r, Ir] = RadialMean(abs(S.amp).^2);
            
            figure, hl = semilogy(r, Ir);
            grid on
            hl.LineWidth = 2;
            xlabel('Radius (pix)')
            ylabel('Mean Intensity')
            hax = gca;
            
            
            
        end % DisplayRadialIntensity
        
        function [Zremap, ampremap, pharemap] = RemapRadial(S, varargin)
            % PIAA radial remapping of amplitude and phase
            % for example, see W:\phaseretrieval\test_remapgsnum_script.m
            % result is saved as amp & phase of S.Eremap
            
            U = CConstants;
            
            % options
            bDebug = CheckOption('debug', false, varargin{:});
            
            hax = [];
            N = size(S.bMask);
            
            % coorindates for this routine are in (mm)
            [x, y, X, Y, R, T] = CreateGrid(N, S.RemapRadialRmax./S.RemapRadialRpix);
          
            % force circular bMask
            % bMaskRemap cannot be larger than bMask
            S.bMaskRemap = (R <= S.RemapRadialRmax) & S.bMask; % (mm)
            
            if bDebug,
                figure_mxn(1,2)
                
                rgbM = zeros([size(S.bMask) 3]);
                rgbM(:,:,1) = 128*S.bMask;
                rgbM(:,:,3) = 128*S.bMaskRemap;
                subplot(1,2,1), imshow(rgbM), title('bMask and bMask for remap')
                
                subplot(1,2,2), plot(X(S.bMaskRemap)/U.MM, Y(S.bMaskRemap)/U.MM, '.')
                axis image
                xlabel('Camera X (mm)'), ylabel('Camera Y (mm)')
                title('Pupil Mask Physical Size at PIAA')
            end
            
            % piaa out map to piaa in
            Rm = zeros(size(R));
            Rm(S.bMaskRemap) = interp1(S.RemapRadialR2, S.RemapRadialR1, R(S.bMaskRemap), 'pchip', nan);
            Xm = Rm.*cos(T); Ym = Rm.*sin(T);
            
            if bDebug,
                figure, plot(Xm(S.bMaskRemap)/U.MM, Ym(S.bMaskRemap)/U.MM, '.'), grid
                axis image
                xlabel('Camera X (mm)'), ylabel('Camera Y (mm)')
                title('PIAA Out to In Map')
            end
            
            % % interpolate phase using complex numbers at xm, ym
            % Ei_r = interp2(X, Y, real(exp(1j*S.phw_ptt)), Xm, Ym);
            % Ei_i = interp2(X, Y, imag(exp(1j*S.phw_ptt)), Xm, Ym);
            % Ei = Ei_r + 1j*Ei_i;
            % pharemap = angle(Ei);
            
            % interpolate unwrapped phase
            pharemap = zeros(size(S.phunwrap));
            pharemap(S.bMaskRemap) = interp2(X, Y, S.phunwrap, Xm(S.bMaskRemap), Ym(S.bMaskRemap));
            
            % mask off any nan's that result from interpolation
            % this shouldn't be necessary, but sometimes it is, need to
            % look into it.
            S.bMaskRemap(isnan(pharemap)) = 0;

            % zernike fit to phase, coordinates are in mm
            Zremap = zernikefit(X(S.bMaskRemap), Y(S.bMaskRemap), pharemap(S.bMaskRemap), 1:11, S.RemapRadialRmax, 'Noll');

            % interpolate amplitude onto remap grid
            ampremap = interp2(X, Y, S.amp, Xm, Ym);
            
            % results S.Eramap, S.phunwrapremap, and S.bMaskRemap are stored in the object
            % properties because CGS is < handle
            S.phunwrapremap = pharemap;
            S.Eremap = ampremap .* exp(1j*pharemap);

        end % RemapRadial
        
        function [hax] = DisplayRemap(S, varargin)
            
            if isempty(S.Eremap),
                RemapRadial(S, varargin{:});
            end
            
            ampremap = S.bMaskRemap .* abs(S.Eremap);
            pharemap = S.bMaskRemap .* S.phunwrapremap;
            %pharemap = S.bMaskRemap .* angle(S.Eremap);
            
            hfig = figure_mxn(2,2);
            
            hax(1,1) = subplot(2,2,1);
            imageschcit(S.amp), colorbartitle('Amplitude'), title(['gsnum ' num2str(S.gsnum) '; Before Remapping'])
            
            hax(1,2) = subplot(2,2,2);
            imageschcit(S.phw_ptt), colorbartitle('Phase (rad)'), title(['gsnum ' num2str(S.gsnum) ';Before Remapping'])
            
            hax(2,1) = subplot(2,2,3);
            imageschcit(ampremap), colorbartitle('Amplitude'), title(['gsnum ' num2str(S.gsnum) ';After Remapping'])
            
            hax(2,2) = subplot(2,2,4);
            imageschcit(pharemap), colorbartitle('Phase (rad)'), title(['gsnum ' num2str(S.gsnum) ';After Remapping'])
            

        end % DisplayRemap
    
        function [ZZ, nzoutnzin, phfit, phresidual, A] = ZernZernRemapFit(S, varargin)            
            % [ZZ, nzoutnzin, phfit, phresidual] = ZernZernRemapFit(CGSobject, varargin)
            %
            % 2020-03-28
            % copied from
            % /home/dmarx/PIAA/hcim_testbed_run100/phaseretrieval_analysis/ZernZernRemapFit.m
            %
            % options:
            % bDisplay = CheckOption('display', false, varargin{:});
            % S.RemapRadialRpix = CheckOption('RemapRadialRpix', S.RemapRadialRpix, varargin{:});
            % bDebug = CheckOption('debug', false, varargin{:});
            % nzout = CheckOption('nzout', 1:4, varargin{:});
            % nzin = CheckOption('nzin', 2:4, varargin{:});
            % poly_order = CheckOption('polyorder', 'Noll', varargin{:});
            
            U = CConstants;
            
            % options
            bDisplay = CheckOption('display', true, varargin{:});
            S.RemapRadialRpix = CheckOption('RemapRadialRpix', S.RemapRadialRpix, varargin{:});
            bDebug = CheckOption('debug', false, varargin{:});
            nzout = CheckOption('nzout', 1:4, varargin{:});
            nzin = CheckOption('nzin', 2:4, varargin{:});
            poly_order = CheckOption('polyorder', 'Noll', varargin{:});
            titlestr = CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});
            
            if isempty(S.Eremap),
                S.RemapRadial;
            end
            
            % % from CGS, apply Zernikes to unwrapped (but not re-mapped)
            % phase, S.phunwrap
            % 
            % piaa out map to piaa in
            % % coorindates for this routine are in (mm)
            [x, y, X, Y, R, T] = CreateGrid(size(S.bMask), S.RemapRadialRmax./S.RemapRadialRpix);
            Rm = zeros(size(R));
            
            % maps back-end to front end
            %Rm(S.bMaskRemap) = interp1(S.RemapRadialR2, S.RemapRadialR1, R(S.bMaskRemap), 'pchip', nan);
            % maps front end to back end
            Rm(S.bMaskRemap) = interp1(S.RemapRadialR1, S.RemapRadialR2, R(S.bMaskRemap), 'pchip', nan);
            Tm = T; % theta is same for PIAA in and out
            
            % matrix equation:
            % phase(:) = [Aout Ain] * [Zout; Zin]
            % phase(r,t) = [P(rout,t) P(rin, t)] * [Zout; Zin]
            % N x 1       N x nzout + Nzin   (Nzout+Nzin) x 1
            P = zernikepolynomials(poly_order);
            Rz = S.RemapRadialRmax;
            M = length(R(S.bMaskRemap));
            A = zeros(M, length(nzout)+length(nzin));
            
            Rp = R(S.bMaskRemap)./Rz; % normalized PIAA out radial coordinates
            Tp = T(S.bMaskRemap);
            for ii = 1:length(nzout),
                A(:,ii) = P{nzout(ii)}( Rp, Tp );
            end
            Rmp = Rm(S.bMaskRemap)./Rz;
            Tmp = Tm(S.bMaskRemap);
            for ii = 1:length(nzin),
                A(:,length(nzout)+ii) = P{nzin(ii)}( Rmp, Tmp );
            end
            
            % take a look at the modes being fitted:
            if bDebug,
                Ntmp = max([length(nzout) length(nzin)]);
                figure_mxn(2, Ntmp)
                for ii = 1:length(nzout)
                    subplot(2, Ntmp, ii)
                    imA = zeros(size(S.bMaskRemap));
                    imA(S.bMaskRemap) = A(:,ii);
                    imageschcit(x/U.MM, y/U.MM, imA)
                    title(['Zout ' num2str(nzout(ii))])
                end
                for ii = 1:length(nzin)
                    subplot(2, Ntmp, Ntmp+ii)
                    imA = zeros(size(S.bMaskRemap));
                    imA(S.bMaskRemap) = A(:, length(nzout)+ii);
                    imageschcit(x/U.MM, y/U.MM, imA)
                    title(['Zin ' num2str(nzin(ii))])
                end
                
                
            end
            
            
            % % L1 norm solution:
            % rhs = S.phunwrap(S.bMaskRemap);
            % f = @(zz) sum(abs(A*zz - rhs));
            % opts = optimset(@fminsearch);
            % opts.MaxFunEvals = 100*opts.MaxFunEvals;
            % ZZL1 = fminsearch(f, zeros(11,1))
            
            % LSE solution:
            if any(isnan(S.phunwrap(S.bMaskRemap))),
                error('check mask regions');
            end
            ZZ = A \ S.phunwrap(S.bMaskRemap);
            
            % best fit phase
            phfit = zeros(size(S.bMaskRemap));
            phfit(S.bMaskRemap) = A*ZZ;
            
            % residual error:
            phresidual = zeros(size(S.bMaskRemap));
            phresidual(S.bMaskRemap) = S.phunwrap(S.bMaskRemap) - A*ZZ;
            
            % return:
            nzoutnzin = [nzout nzin];
            
            %%%%%% plot results
            if bDisplay,
                figure_mxn(2,3)

                subplot(2,3,1), imageschcit(S.x, S.y, abs(S.E))
                title([titlestr '; Amplitude'])
                
                %subplot(2,3,2), imageschcit(S.phw_ptt.*S.bMaskRemap)
                subplot(2,3,2), imageschcit(S.x, S.y, mod2pi(S.phunwrap).*S.bMaskRemap)
                colorbartitle('Phase (rad)')
                set(gca,'clim',pi*[-1 1])
                title([titlestr '; Phase rms\phi = ' num2str(S.rmsPha,'%.3f')])               
                
                subplot(2,3,3), imageschcit(S.x, S.y, phresidual), 
                colorbartitle('Phase (rad)')
                rmse = rms(phresidual(S.bMaskRemap));
                title(['Residual Fit, rms error = ' num2str(rmse,'%.2f') 'rad'])
                
                % bar graph zernike order, don't plot piston
                %figure,
                subplot(2,3,4), imageschcit(S.x, S.y, phfit), 
                colorbartitle('Phase (rad)'), title('Fit Phase')
                
                
                subplot(2,3,5:6)
                hh = bar(ZZ(2:end)); grid
                set(gca,'XTick', 1:length([nzout(2:end) nzin]))
                set(gca, 'XTickLabel', num2str([nzout(2:end) nzin]'))
                %set(gca,'Position',[0.05 0.15 0.9 0.8])
                ylabel('Zernike Coeff (rms rad)')
                xlabel('Zernike # [PIAA Out][PIAA In] (Noll Order)')
                % set(gcf,'Position', [0 0 1600 400] + [1 1 0 0].*get(gcf,'Position'))
                % title(['gsnum ' num2str(gsnum) ]); %'; Foc Ref ' sz])
                set(gca,'ylim',2*[-1 1])
                % set(gcf,'Position',0.6*get(gcf,'Position'))
                
            end
            

        end % ZernZernRemapFit
    
        function S = AdjustUnrapRegionPiston(S)
            
            % remove discontinuites across region boundaries if pupil mask has separate regions
            % find objects
            % eliminate isolated pixels, etc

            [B,L,N,A] = bwboundaries(S.bMask, 'noholes');

            bwmask = S.bMask;
            for il = 1:N,
                objarea(il) = bwarea(L==il);
                if objarea(il) < 100,
                    %S.bMask(L==il) = true;
                    bwmask(L==il) = 0;
                end
            end
            [B,L,N,A] = bwboundaries(bwmask, 'noholes');
            figure, imageschcit(1,1,L), hold on, plot(B{1}(:,2), B{1}(:,1), '-r')
            
            %
            [x, y, X, Y] = CreateGrid(bwmask, 1, 1, 'origin', '1-offset');
            
            % find approximate inner and outer diameters
            maskcom = calcCOM(x, y, bwmask);
            [Xc, Yc] = meshgrid(x-maskcom(1),y-maskcom(2)); R = hypot(Xc, Yc);
            figure, HR = histogram(R(bwmask));
            
            radius_id = HR.BinEdges(1)+0.5*HR.BinWidth;
            radius_od = HR.BinEdges(end-1)+0.5*HR.BinWidth;
            rpad = 0.1*(radius_od - radius_id);
            bwmask(R>(radius_od-rpad)) = 1;
            bwmask(R<(radius_id+rpad)) = 1;
            bwmask = ~bwmask;
            [Bs, Ls, Ns, As] = bwboundaries(bwmask, 'noholes'); % struts
            figure, imageschcit(1,1,Ls), colorbar, hold on, plot(Bs{1}(:,2), Bs{1}(:,1), '-r')
            if Ns ~= N, error([num2str(N) ' mask regions, but ' num2str(Ns) ' struts']); end
            
            % for each struct, find neighboring mask regions and check phase delta across the strut
            for ns = 1:Ns,
                for nn = 1:N,
                    dmin_ii = zeros(length(Bs{ns}(:,1)),1);
                    for ii = 1:length(Bs{ns}(:,1)),
                        dmin_ii(ii) = min(hypot(Bs{ns}(ii,1)-B{nn}(:,1), Bs{ns}(ii,2)-B{nn}(:,2)));
                    end
                    dmin_nn(nn) = min(dmin_ii);
                end % for each pupil region
                % there should be two regions with dmin close to 1 pixel
                [dmin, nnmin] = sortdata({dmin_nn, 1:N});
                if dmin(2) > 2, error('something is wrong'); end
                
                % then neighboring mask regions are nnmin(1:2)
                % in clockwise order:
                nnmin = sort(nnmin(1:2));
                
                % mean phase of each neighboring region near the boundary
                % center of this strut (not centroid, just simple center of the boundary points):
                strut_xc = mean(Bs{ns}(:,2));
                strut_yc = mean(Bs{ns}(:,1));
                
                % find pixels in each region near the strut
                for inn = 1:2,
                    [ddd, xxx, yyy] = sortdata({hypot(X(L==nnmin(inn))-strut_xc, Y(L==nnmin(inn))-strut_yc), X(L==nnmin(inn)), Y(L==nnmin(inn))});
                    iuse = 1:ceil(0.05*length(ddd));
                    meanph(inn) = mean(S.phunwrap(sub2ind(size(S.phunwrap), yyy(iuse),xxx(iuse))));
                    hold on, if inn==1, plot(xxx(iuse),yyy(iuse),'or'), else, plot(xxx(iuse),yyy(iuse),'xr'), end
                end
                
                deltaph = abs(diff(meanph));
                if deltaph > pi,
                    % adjust the region farther from zero towards zero
                    if abs(meanph(1)) > abs(meanph(2)),
                        inn = 1;
                    else
                        inn = 2;
                    end
                    
                    psign = sign(meanph(inn));
                    S.phunwrap(L==nnmin(inn)) = S.phunwrap(L==nnmin(inn)) - psign*ceil((deltaph-pi)/(2*pi))*2*pi;
                    
                end % if deltaph > pi
                
            end % for each strut
            
    
        end % AdjustUnwrapRegionPiston
        
    end % methods
    
end % classdef


function phw_ptt = RemovePTTZ(phw, bMask)
[x, y, X, Y, R] = CreateGrid(phw);
Rz = max(R(bMask));
Z = zernikefit(X(bMask), Y(bMask), phw(bMask), 3, Rz);
phw_ptt = zeros(size(X));
phw_ptt(bMask) = phw(bMask) - zernikeval(Z, X(bMask), Y(bMask), Rz);
%phw_ptt(:) = mod2pi(phw(:) - zernikeval(Z, X(:), Y(:), Rz));

end % RemovePTTZ