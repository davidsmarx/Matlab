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
    %                     case 'mcb_alllens'
    %                     case 'mcb_twolens'
    %                     case 'piaacmc'
    %                     case 'omc_mswc'
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
        nRot90     % rotate all fits images after reading using rot90(Im, nRot90) (default = no rotate)
        zAmpPlanes % camz for each amp image
        zunits = 1;% units for reading zAmpPlanes from fits
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

        params % from yaml
        
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

            wavelength_kwd = 'lam'; % default value
            wavelength_units = 1; % default = m

            % default value for bn
            if ~exist('bn','var') || isempty(bn),
                %bn = '/home/dmarx/HCIT/DST/phaseretrieval_20180605/reduced/gsdst_';
                %bn = '/proj/dst/data/dB_PR/gsdst_';
                %bn = 'omc_mswc';
                bn = 'cgi_fft';
            end
            
            switch lower(bn),
                case 'dst'
                    bn = ['/proj/dst/data/hcim/dB_PR/gsdst_' num2str(gsnum,'%d')];
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
                    bn = ['/home/dmarx/HCIT/SPC_disc/gsspc_20171204/reduced/gsspc_' num2str(gsnum)];
                    
                case 'mcb_spc'
                    bn = ['/home/dmarx/HCIT/MCB_SPC/phaseretrieval/reduced/gsomc_no' num2str(gsnum,'%04d')];
                    
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
                    wavelength_units = U.MM;
                    
                    %
                    S.zunits = U.MM;
                    
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
                    %if ~isempty(varargin),  trialname = [varargin{1} '_']; end
                    if ~isempty(varargin),  trialname = [varargin{1}]; end
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
                    
                    wavelength_kwd = 'lam';
                    
                case 'omc_mswc'
                    trialname = CheckOption('trialname', '', varargin{:});
                    
                    bn = ['/home/hcit/OMC/phaseretrieval/reduced/prout_' trialname num2str(gsnum,'%03d')];
                    
                    % get dir listing of raw camera images
                    S.listPupImDir = dir(PathTranslator(...
                        ['/proj/mcb/data/excam/*/pr_gs_' num2str(gsnum,'%04d') '/emccd.*.fits']...
                        ));
                    S.listSrcImDir = dir(PathTranslator(...
                        ['/proj/mcb/data/excam/*/pr_par_' num2str(gsnum,'%04d') '/emccd.*.fits']...
                        ));
                    
                    wavelength_kwd = 'lam';
                    
                case 'cgi_fft'
                    
                    % Note: add check to see if data already exists
                    % locally, then skip scp

                    % scp reduced files from yzma to local
                    % first make local folder
                    base_name = ['prnum_' num2str(gsnum, '%06d')];
                    local_path = PathTranslator('~/WFIRST/VA_FFT_activities/FFT/pr/FFT_reduced/');
                    if ~exist(fullfile(local_path, base_name), 'dir')
                        % get the pr data package from kronk
                        url = ['https://kronk.jpl.nasa.gov:8000/flight/pr/' num2str(gsnum, '%06d') '.zip'];
                        zip_fn = websave(fullfile(local_path, [base_name '.zip']), url);
                        list_fn = unzip(zip_fn, local_path);
                    end
                    % else
                    % apparently this reduced data already transferred
                    
                    bn = fullfile(local_path, base_name, filesep);
                    
                    wavelength_kwd = 'lam';

                otherwise
                    % do nothing, bn is explicit, check that it is
                    % valid path
                    bn = [bn num2str(gsnum,'%d')];
                    %disp(bn);
                    if isempty(ls(PathTranslator([bn '*'])))
                        error(['invalid bn: ' bn]);
                    end
                    
            end % switch bn
            
            %disp(['opening: ' PathTranslator([bn num2str(gsnum,'%03d') 'amp.fits'])]);
            try
                %ampinfo = fitsinfo(PathTranslator([bn num2str(gsnum,'%03d') 'amp.fits']));
                ampinfo = fitsinfo(PathTranslator([bn 'amp.fits']));
            catch ME
                %disp(ME.message);
                fprintf('failed to open fits file: \n%s\n',PathTranslator([bn 'amp.fits']));
            end
            
            % check options
            S.nRot90 = CheckOption('rot90', 0, varargin{:});
 
            S.gsnum = gsnum;
            S.bn = bn;
            S.amp = rot90( fitsread(PathTranslator([bn 'amp.fits'])), S.nRot90);
            S.ph = rot90( fitsread(PathTranslator([bn 'ph.fits'])), S.nRot90);  % unwrapped phase
            S.phw = rot90( fitsread(PathTranslator([bn 'phwrap.fits'])), S.nRot90); % = angle(eref), wrapped phase
            S.amp_keys = ampinfo.PrimaryData.Keywords;
            S.wavelength = FitsGetKeywordVal(S.amp_keys, wavelength_kwd)*wavelength_units;

            % real parms yaml file
            try
                list_fn = {PathTranslator([S.bn 'parms.yml']), PathTranslator([S.bn 'parms.yaml'])};
                if ~any(isfile(list_fn))
                    error('No parameter yaml file found');
                end
                for ifn = 1:length(list_fn)
                    if isfile(list_fn{ifn})
                        S.params = yaml.loadFile(list_fn{ifn});
                        break;
                    end
                end
                
            catch ME
                disp(ME.message);
            end
            
            % phase unwrap in WFSC puts the mask used for PR in the second
            % hdu (starting Feb 2021)
            finfo = fitsinfo(PathTranslator([bn 'ph.fits']));
            if length(finfo.Contents) > 1,
                fn = [bn 'ph.fits'];
                S.phunwrap_bMask = logical(rot90( fitsread(PathTranslator(fn),'image'), S.nRot90) );
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
            S.RemovePTTfft; % creates first estimate of S.phw_ptt using FFT
            %             S.phunwrap = RemovePTTZ(S.phunwrap, S.bMask);
            %             S.phw_ptt  = mod2pi(S.phunwrap);

            % unwrap phase using better unwrap routine, but requires mask
            %phw = S.phw;
            phw = S.phw_ptt + 1 - 1;
            phw(~S.bMask) = NaN;
            S.phunwrap = unwrap_phase(phw);
            S.phunwrap(~S.bMask) = 0;
            S.phw_ptt = RemovePTTZ(S.phw_ptt, S.bMask); % fine-tune using Zernike

            S.phunwrap = S.phw_ptt;

            % S.E
            %S.E = S.amp .* exp(1i*S.phw_ptt);
            S.E = S.amp .* exp(1i*S.phw);
        
            % remove tip/tilt from the given phase unwrapped
            [~, phaimg, ~, ~] = ZernikeAnalysis(S.ph, 'modes', 1:3, 'bMask', S.bMask, 'isphase', true);
            S.ph = phaimg;
            
            % load the radial mapping
            % should be part of the bn switch
            % see email from Dan Sirbu "RE double check r2 r1 definition"
            if isequal(bn, 'piaacmc')
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
            
            S.cAmpPlanes{1} = rot90( fitsread(PathTranslator(cmp_fn)), S.nRot90);
            S.zAmpPlanes(1) = FitsGetKeywordVal(finfo.PrimaryData.Keywords, 'Z')*S.zunits;
            for ii = 2:NIm,
                S.cAmpPlanes{ii} = rot90( fitsread(PathTranslator(cmp_fn),'image',ii-1), S.nRot90 );
                S.zAmpPlanes(ii) = FitsGetKeywordVal(finfo.Image(ii-1).Keywords, 'Z')*S.zunits;
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
            % phplot = CheckOption('phplot', 'phw_ptt', varargin{:}); %
            %     other choices = 'phw_ptt', 'ph', 'phw', 'phunwrap', 'angleE'
            % ampplot = CheckOption('ampplot', 'absE', varargin{:});
            %     other choices 'amp'
            % stitle = CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});
            % bRemoveTipTilt = CheckOption('removetiptilt', true, varargin{:});            
            % dph_units = CheckOption('dph_units', 1, varargin{:}); % default = 'rad' (radians), 'nm', 'waves', or double
            % dph_units_str = CheckOption('dph_units_str', 'Phase (rad)', varargin{:});

            U = CConstants;
            
            pMask = CheckOption('pMask', S.bMask, varargin{:});
            xylim = CheckOption('xylim', 1.1*max(S.R(S.bMask)), varargin{:});
            climph = CheckOption('climph', [], varargin{:});
            phplot = CheckOption('phplot', 'ph', varargin{:}); % or S.(phplot), e.g. 'angleE'
            ampplot = CheckOption('ampplot', 'absE', varargin{:});
            stitle = CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});
            bRemoveTipTilt = CheckOption('removetiptilt', true, varargin{:});
            dph_units = CheckOption('dph_units', 1, varargin{:}); % default = 'rad' (radians), 'nm', 'waves', or double
            dph_units_str = CheckOption('dph_units_str', 'rad', varargin{:});
            
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

            % remove tip/tilt?

            % bMask is only for phase plot
            % check options for what to plot
            % [ZZ, rz, pharesidual] = ZernikeFit(S, nz, varargin)
            %
            % zernike fit using bMask pixels
            % nz = array of zernike modes to fit
            %
            % phase = CheckOption('phase', 'ph', varargin{:});
            % bDisplay = CheckOption('display', true, varargin{:});
            %[~, ~, ph] = S.ZernikeFit(1:3, 'phase', phplot, 'display', false);
            
            switch phplot
                case 'angleE'
                    ph = angle(S.E);
                otherwise
                    ph = S.(phplot);
            end
            

            if ischar(dph_units),
                switch dph_units
                    case 'nm'
                        if ~isempty(S.wavelength),
                            dph_units = U.NM./(S.wavelength./(2*pi));
                            dph_units_str = 'nm';
                        else
                            warning('wavelength empty, dph units = radians');
                            dph_units = 1;
                        end
                    case {'wave', 'waves'},
                        dph_units = 2*pi;
                        dph_units_str = 'waves';
                        
                    case 'rad'
                        dph_units = 1;
                        dph_units_str = 'rad';
                    otherwise
                        error(['unknown phase units: ' dph_units]);
                end
            end % ischar(dph_units)
             

            if isempty(pMask), pMask = ones(size(S.E)); end
            hax(2) = subplot(1,2,2);
            imageschcit(S.x, S.y, pMask.*ph/dph_units)
            colorbartitle(['WFE (' dph_units_str ')'])
            set(gca,'xlim',xylim*[-1 1],'ylim',xylim*[-1 1])
            if isempty(climph),
                climph = AutoClim(ph/dph_units, 'symmetric', true);
            end
            set(gca,'clim',climph)
            title(['rms\phi = ' num2str(rms(ph(pMask))/dph_units,'%.3f') ' ' dph_units_str])
            
            
            
        end % DisplayGS

        function [hfig, hax, dphaResult] = DisplayGSrefGS(S, Sref, varargin)
            % [hfig, hax, dphaResult] = DisplayGSrefGS(S, Sref, options)
            %
            % options:
            %    ('hfig', figure_mxn(2,2), varargin{:});
            %    ('usebMask', true, varargin{:});
            %    ('pMask', S.bMask, varargin{:}); % bMask only for phase display
            %    ('removeDefocus', false, varargin{:});
            %    ('removeZ2Z3', true, varargin{:});
            %    ('doRegister', false, varargin{:}); % (false), true
            %    ('phplot', 'ph', (default) 'phw_ptt', 'phunwrap', 'angleE' S.(phplot)
            %    ('xylim', 1.1*max(S.R(S.bMask)), varargin{:});
            %    ('dphclim', [], varargin{:});
            %    ('climph', [], varargin{:});
            %    ('dph_units', 1, varargin{:}); % default = radians, 'nm', 'waves', or double
            %    ('dph_units_str', 'Phase (rad)', varargin{:});
            %   
            % dphaResult = struct(...
            %    'ZZ', ZZ ...  % only when removeDefocus is true
            %    ,'phaimg', phaimg ... % dphase with Zernike 1:3 = 0
            %    ,'dpha', dpha ... % residual dphase with Zernike 1:4 = 0
            %    ,'sOptions', sOptions ... % only when removeDefocus is true
            %    );

            U = CConstants;
            
            % parse options
            hfig = CheckOption('hfig', [], varargin{:});
            usebMask = CheckOption('usebMask', true, varargin{:});
            pMask = CheckOption('pMask', S.bMask, varargin{:}); % bMask only for phase display
            removeDefocus = CheckOption('removeDefocus', false, varargin{:});
            removeZ2Z3 = CheckOption('removeZ2Z3', true, varargin{:});
            doRegister = CheckOption('doRegister', false, varargin{:});
            phplot = CheckOption('phplot', 'ph', varargin{:}); % S.(phplot)
            xylim = CheckOption('xylim', [], varargin{:});
            climdph = CheckOption('dphclim', [], varargin{:});
            climph = CheckOption('climph', [], varargin{:});
            dph_units = CheckOption('dph_units', 1, varargin{:}); % default = radians, 'nm', 'waves', or double
            dph_units_str = CheckOption('dph_units_str', 'Phase (rad)', varargin{:});
            dph_format = CheckOption('dph_format', '%.3f', varargin{:}); % format for title string

            % which phase map to plot
            switch phplot
                case 'angleE'
                    funPhPl = @(S) angle(S.E);
                otherwise
                    funPhPl = @(S) S.(phplot);
            end
            
            if ishandle(hfig),
                figure(hfig);
            else
                hfig = figure_mxn(2,2);
            end

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
                        ,'amp', CropImage(Sref.amp, [], [0 0], nc, nr) ...
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
            title(['gsnum ' num2str(S.gsnum) ' -  Ref ' num2str(Sref.gsnum) ' Amp'])
            
            
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
            elseif removeZ2Z3
                [ZZ, phaimg, dpha, sOptions] = ZernikeAnalysis(funPhPl(S) - Sreftmp.phplot, ...
                    'modes', 1:3, 'bMask', S.bMask, 'isphase', true);
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
                dpha = pMask .* dpha;
            end
            if ischar(dph_units),
                switch dph_units
                    case 'nm'
                        if ~isempty(S.wavelength),
                            dph_units = U.NM./(S.wavelength./(2*pi));
                            dph_units_str = 'WFE (nm)';
                            dph_format = '%.1f';
                        else
                            warning('wavelength empty, dph units = radians');
                            dph_units = 1;
                        end
                    case {'wave', 'waves'},
                        dph_units = 2*pi;
                        dph_units_str = 'WFE (waves)';                        
                        
                    otherwise
                        error(['unknown phase units: ' dph_units]);
                end
            end % ischar(dph_units)
            
            him = imageschcit(S.x, S.y, dpha./dph_units);
            colorbartitle(dph_units_str)
            set(gca,'xlim',xylim*[-1 1],'ylim',xylim*[-1 1])
            if isempty(climdph), set(gca,'clim',AutoClim(dpha./dph_units,'symmetric',true,'pctscale',100))
            else set(gca,'clim',climdph)
            end
            title(['gsnum ' num2str(S.gsnum) ' Ref ' num2str(Sref.gsnum) ', rms \Delta = ' num2str(rms(dpha(pMask))/dph_units,dph_format) dph_units_str])
            
            % append more to dphaResult
            dphaResult.dpharms = rms(dpha(pMask));
            dphaResult.pMask = pMask;
            
        end % DisplayGSrefGS

        function phw_ptt = RemovePTTfft(S)
            % remove large scale piston, tip, tilt by translating the fft
            % s.t. the fft peak is at the center
            
            Zs = S.amp.*exp(1i*S.phw);

            %ZZ = fftshift(fft2(fftshift(S.bMask.*Zs)));
            ZZ = fftshift(fft2(fftshift(Zs)));
            
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
        
        function [ZZout, rz, pharesidual, sFitParms] = ZernikeFit(S, nz, varargin)
            % [ZZ, rz, pharesidual, sFitParms] = ZernikeFit(S, nz, varargin)
            %
            % zernike fit using bMask pixels
            % nz = array of zernike modes to fit (default = 1:11)
            %
            % phase = CheckOption('phase', 'ph', varargin{:});
            % bmask = CheckOption('bmask', S.bMask, varargin{:});
            % bDisplay = CheckOption('display', true, varargin{:});
            % titlestr = CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});
            % xylim = CheckOption('xylim', [], varargin{:});
            % phresclim = CheckOption('phresclim', [], varargin{:});
            % ylimZ = CheckOption('ylimzplot', [], varargin{:});
            % CheckOption('zernikeunits', 'nm', varargin{:}); % 'nm', 'rad'
            %
            % return:
            % ZZout = zernike coefficients, ZZout(1:3) are always piston,
            %     tip, tilt
            % rz = normalization radius (pixels)
            % pharesidual = S.(phase) - zernikeval(ZZout)

            % default nz
            if ~exist('nz','var')
                nz = 1:11;
            end
            
            phasefieldname = CheckOption('phase', 'ph', varargin{:});
            bmask = CheckOption('bmask', S.bMask, varargin{:});
            bDisplay = CheckOption('display', true, varargin{:});
            hfig = CheckOption('hfig', [], varargin{:});
            titlestr = CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});
            xylim = CheckOption('xylim', [], varargin{:});
            phresclim = CheckOption('phresclim', [], varargin{:});
            ylimZ = CheckOption('ylimzplot', [], varargin{:});
            zzunits = CheckOption('zernikeunits', 'nm', varargin{:}); % 'nm' or 'rad'
            
            % fit should always include piston, tip, tilt, even if not
            % included in requested nz
            nzfit = nz(:);
            if ~any(nzfit == 3), nzfit = [3; nzfit]; end
            if ~any(nzfit == 2), nzfit = [2; nzfit]; end
            if ~any(nzfit == 1), nzfit = [1; nzfit]; end
            
            %
            rz = max(S.R(bmask));
            %             ZZ = zernikefit(S.X(bmask), S.Y(bmask), S.(phasefieldname)(bmask), nzfit, rz, 'noll');
            %             phwfit = nan(size(S.X));
            %             phwfit(bmask) = zernikeval(ZZ, S.X(bmask), S.Y(bmask), rz, 'noll', 'nz', nzfit);
            %             pharesidual = nan(size(S.X));
            %             pharesidual(bmask) = mod2pi(S.(phasefieldname)(bmask) - phwfit(bmask));
            [ZZ, phaimg_1_3, pharesidual, sFitParms] = ZernikeAnalysis(S.(phasefieldname),...
                'isphase', true, 'bMask', bmask, 'Rnorm', rz, 'modes', nzfit, 'polyorder', 'Noll',...
                'do_phaseunwrap', false);
            phwfit = bmask.*sFitParms.phafit_ptt;
            
            % return zernike coeffs in requested units
            switch zzunits
                case 'rad'
                    % no change to ZZ
                    ZZout = ZZ;
                    strUnits = '(rms rad)';
                
                case 'nm'
                    ZZout = S.wavelength*ZZ/(2*pi)/1e-9;
                    strUnits = '(rms nm)';                    
                    
                otherwise
                    warning(['invalid units: ' zzunits ', default to radians']);
                    strUnits = '(rms rad)';    
            end
            
            %%%%%% plot results
            if bDisplay,

                % determine plot width
                if isempty(xylim),
                    xylim = 1.1*max(S.R(bmask));
                    xylim = 5*ceil(xylim/5.0); % to the nearest multiple of 5
                    xylim = xylim*[-1 1];
                end

                if isempty(hfig),
                    hfig = figure_mxn(2,3);
                else
                    figure(hfig)
                end

                subplot(2,3,1), imageschcit(S.x, S.y, abs(S.E))
                title([titlestr '; Amplitude'])
                set(gca,'xlim',xylim,'ylim',xylim);

                % phase, remove 1 to 3
                subplot(2,3,2), imageschcit(S.x, S.y, phaimg_1_3.*bmask) % mod2pi(S.phunwrap).*bmask)
                colorbartitle('Phase (rad)')
                %set(gca,'clim',pi*[-1 1])
                title([titlestr '; Phase rms\phi = ' num2str(rms(phaimg_1_3(bmask)),'%.3f')])               
                set(gca,'xlim',xylim,'ylim',xylim);
                
                % residual phase
                subplot(2,3,3), imageschcit(S.x, S.y, pharesidual),
                colorbartitle('Phase (rad)')
                rmse = rms(pharesidual(bmask));
                title(['Residual Fit, rms error = ' num2str(rmse,'%.3f') 'rad'])
                set(gca,'xlim',xylim,'ylim',xylim);
                if ~isempty(phresclim), set(gca,'clim',phresclim); end
                
                % fit phase
                subplot(2,3,4), imageschcit(S.x, S.y, phwfit), 
                colorbartitle('Phase (rad)'), title('Fit Phase')
                set(gca,'xlim',xylim,'ylim',xylim);                
                
                % bar graph zernike order, don't plot piston
                subplot(2,3,5:6)
                [Zplot, NZplot] = filterdata(nzfit > 3, ZZout, nzfit);
                hh = bar(NZplot, Zplot); grid
                ylabel(['Zernike Coeff ' strUnits])
                xlabel('Zernike # (Noll Order)')

                if ~isempty(ylimZ),
                    set(gca,'ylim',ylimZ)
                end
                
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
            
            % weights included?
            [~, ~, nw] = size(S.cAmpPlanes{1});
            
            CCor = @(a,b) a(:)'*b(:)./sqrt( (a(:)'*a(:)) * (b(:)'*b(:)) );
            CCorq = @(q) CCor(q(:,:,1),q(:,:,2));
            CCint = @(q) CCor(abs(q(:,:,1)).^2, abs(q(:,:,2)).^2);
            %CCwint= @(q) CCor(q(:,:,4).*abs(q(:,:,1)).^2, q(:,:,4).*abs(q(:,:,2)).^2);
            CCwint = @(q) S.WeightedIntensityCorrelation(q(:,:,4), abs(q(:,:,1)).^2, abs(q(:,:,2)).^2);
            
            NIm = length(S.cAmpPlanes);
            if nw == 4,
                cc = zeros(NIm,3);
            else
                cc = zeros(NIm,2);
            end

            for ii = 1:NIm,
                cc(ii,1) = CCorq(S.cAmpPlanes{ii});
                cc(ii,2) = CCint(S.cAmpPlanes{ii});
                if nw == 4,
                    cc(ii,3) = CCwint(S.cAmpPlanes{ii});
                end
            end

            % this is the merit function used to parameter searching
            MF = [length(cc(:,1)) - sum(cc(:,1)) length(cc(:,2))-sum(cc(:,2))];
            if nw == 4,
                MF = [MF length(cc(:,3))-sum(cc(:,3))];
            end
            
        end % AmpCorrMetric
        
        function cc = WeightedIntensityCorrelation(S, ww, intensity_1, intensity_2)
            CCor = @(a,b) a(:)'*b(:)./sqrt( (a(:)'*a(:)) * (b(:)'*b(:)) );

            ww(isnan(ww)) = 0;
            cc = CCor(ww.*intensity_1, ww.*intensity_2);
        end
        
        function [wsqrt, imvars] = CalcWeights(S)
            %
            %
            
            if isempty(S.cAmpPlanes),
                S.ReadAmpImages;
            end
            
            rn = 2.6;
            gaine = 0.5;
            
            for iim = 1:length(S.cAmpPlanes),
                Iimg = squeeze(S.cAmpPlanes{iim}(:,:,1)).^2; % measured intensity
                imvars{iim} = (rn.^2 + Iimg/gaine)./(sum(Iimg(:)).^2);
                
                amp_e = squeeze(S.cAmpPlanes{iim}(:,:,2)); % estimated amplitude
                alpha_2 = sum(amp_e(:).^2);
                
                wsqrt{iim} = sqrt(1./(alpha_2.*imvars{iim}));
                
            end
            
        end
        
        function [hfig, hax] = DisplayAmpPlane(S, ipl, varargin)
            % [hfig, hax] = DisplayAmpPlane(S, list_ipl)
            %
            % compare camera and estimated
            % options:
            %   Option('value', 'amp' (default), 'intensity'
            %   Option('clim', [], varargin{:});
            %   Option('scale', 'log', varargin{:}); % (log), linear
            %   Option('xylim', [], varargin{:});
            %   Option('type', 'tile', varargin{:}); % or 'blink' for ImageCube to blink
            %   Option('hfig', [], varargin{:}); %
            
            if isempty(S.cAmpPlanes), S.ReadAmpImages; end
            
            plAmpOrInt = CheckOption('value', 'amp', varargin{:});
            clim = CheckOption('clim', [], varargin{:});
            scale = CheckOption('scale', 'log', varargin{:});
            xylim = CheckOption('xylim', [], varargin{:});
            display_type = CheckOption('type', 'tile', varargin{:});
            hfig = CheckOption('hfig', [], varargin{:});
            
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

            switch display_type
                case 'blink'
                    Nplanes = length(ipl);
                    
                    %imcube = zeros([2*Nplanes, size]);
                    if isempty(hfig), hfig = figure; else, figure(hfig); end
                    
                    for iipl = 1:Nplanes,
                        imcube(2*iipl-1,:,:) = funPlot(S.cAmpPlanes{ipl(iipl)},1);
                        imcube(2*iipl,:,:)   = funPlot(S.cAmpPlanes{ipl(iipl)},2);
                    end
                    
                    [hfig, hax] = ImageCube(imcube, 1:2*Nplanes); % 'fTitleStr', @(isl) {['Camera, ipl = ' num2str(ipl(ceil(isl/2)))]
                    %set(hax,'clim',get(gca,'clim'))
                    if ~isempty(xylim)
                        set(hax, 'xlim', xylim, 'ylim', xylim)
                    end
                    
                case 'tile',
                    Nplanes = length(ipl);
                    
                    % create figure if necessary
                    if isempty(hfig), hfig = figure_mxn(2, Nplanes); else, figure(hfig); end
                    
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
                otherwise
                    error(['unknown display type: ', display_type]);
                    
            end % switch display type
            
        end % DisplayAmpPlane
        
        function [hax, Imgs, ha] = DisplayAllPlanes(S, varargin)
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
                    Imlog = logImage(Im);
                    imageschcit(Imlog)
                    %set(gca,'clim',[-4 0] + max(Imlog(:)))
                    
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
        
        function [r, Ir, hfig, hax, hl] = DisplayRadialIntensity(S, ipl, varargin)
            % [r, Ir, hax] = DisplayRadialIntensity(S, ipl, varargin)
            % display mean intensity vs radius
            
            if isempty(S.cAmpPlanes), S.ReadAmpImages; end
            
            [xx, yy] = CreateGrid(S.cAmpPlanes{8}(:,:,1));
            [r_cam, Ir_cam] = RadialMean(xx, yy, abs(S.cAmpPlanes{ipl}(:,:,1)).^2);
            [r_est, Ir_est] = RadialMean(xx, yy, abs(S.cAmpPlanes{ipl}(:,:,2)).^2);
            
            figure, hl = semilogy(r_cam, Ir_cam./max(Ir_cam(:)), r_est, Ir_est./max(Ir_est(:)));
            grid on
            set(hl, 'LineWidth', 2);
            xlabel('Radius (pix)')
            ylabel('Mean Intensity')
            legend('Camera', 'Estimated')
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
            % CheckOption('display', false, varargin{:});
            % CheckOption('RemapRadialRpix', S.RemapRadialRpix, varargin{:});
            % CheckOption('debug', false, varargin{:});
            % CheckOption('nzout', 1:4, varargin{:});
            % CheckOption('nzin', 2:4, varargin{:});
            % CheckOption('polyorder', 'Noll', varargin{:});
            % CheckOption('xylim', 1.1*max(S.R(S.bMask)), varargin{:});
            % CheckOption('title', ['gsnum ' num2str(S.gsnum)], varargin{:});
            
            U = CConstants;
            
            % options
            bDisplay = CheckOption('display', true, varargin{:});
            S.RemapRadialRpix = CheckOption('RemapRadialRpix', S.RemapRadialRpix, varargin{:});
            bDebug = CheckOption('debug', false, varargin{:});
            nzout = CheckOption('nzout', 1:4, varargin{:});
            nzin = CheckOption('nzin', 2:4, varargin{:});
            poly_order = CheckOption('polyorder', 'Noll', varargin{:});
            xylim = CheckOption('xylim', 1.1*max(S.R(S.bMask)), varargin{:});
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
                set(gca, 'xlim', xylim*[-1 1], 'ylim', xylim*[-1 1]);
                title([titlestr '; Amplitude'])
                
                %subplot(2,3,2), imageschcit(S.phw_ptt.*S.bMaskRemap)
                subplot(2,3,2), imageschcit(S.x, S.y, mod2pi(S.phunwrap).*S.bMaskRemap)
                set(gca, 'xlim', xylim*[-1 1], 'ylim', xylim*[-1 1]);
                colorbartitle('Phase (rad)')
                set(gca,'clim',pi*[-1 1])
                title([titlestr '; Phase rms\phi = ' num2str(S.rmsPha,'%.3f')])               
                
                subplot(2,3,3), imageschcit(S.x, S.y, phresidual), 
                set(gca, 'xlim', xylim*[-1 1], 'ylim', xylim*[-1 1]);
                colorbartitle('Phase (rad)')
                rmse = rms(phresidual(S.bMaskRemap));
                title(['Residual Fit, rms error = ' num2str(rmse,'%.2f') 'rad'])
                
                %figure,
                subplot(2,3,4), imageschcit(S.x, S.y, phfit), 
                set(gca, 'xlim', xylim*[-1 1], 'ylim', xylim*[-1 1]);
                colorbartitle('Phase (rad)'), title('Fit Phase')
                
                % bar graph zernike order, don't plot piston                
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