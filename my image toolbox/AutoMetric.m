function [sResults, bMask] = AutoMetric(A, B, varargin)
% [sResults, mask] = AutoMetric(A)
% [sResults, mask] = AutoMetric(A, B)
% [sResults, mask] = AutoMetric(A, B, sOptions)
% 
% A and B are complex arrays of the same size
% compare A and B, separate into 'inside' and 'outside' regions
% compare difference between A and B in each region and test difference
% against eps.
%
% [sResults, mask] = AutoMetric(A)
%     mask = AutoThreshold(A)
%     sResults.thresh = threshold used to make mask
%
% default options
%     sOptions = struct(...
%         'Wgt_mask', [] ... % option to supply a real valued array to create the mask
%         ,'image_type', 'pupil' ... % 'pupil' or 'psf'
%         ,'eps_inside', 1.0e-15 ...
%         ,'eps_outside', 1.0e-12 ...
%         ,'AutoThreshold_Nbins', 21 ...
%         ,'bScaleAmp', false ...
%         ,'logPSF', false ...
%         ,'debug', false ...
%         ,'PSF_thresh_nsig', 2 ...
%         );
%
% if sOptions.Wgt_mask is empty:
%    test 1:
%       mask_a = AutoThreshold(A)
%       mask_b = AutoThreshold(B)
%       if mask_a == mask_b then bPassMask = true
%    test 2:
%       mf = max( abs(A(mask)-B(mask)) ./ min(abs(A(mask),abs(B(mask))) )
%       bPassInside  = mf <= sOptions.eps_inside
%       bPassOutside = mf <= sOptions.eps_outside
%
% if sOptions.Wgt is given:
%    test 1: always pass
%    test 2:
%       mask = AutoThreshold(Wgt)
%       mf = max( abs(A(mask)-B(mask)) ./ min(abs(A(mask),abs(B(mask))) )
%       bPassInside  = mf <= sOptions.eps_inside
%       bPassOutside = mf <= sOptions.eps_outside

% version 2016-03-31
%   added AutoThresholdPSF() for PSF like images where there is not a clear
%   separation in the intensity histogram to separate signal from
%   background
% version 2017-07-19
%   return thresh = threshold used to make mask in all cases.
% version 2017-08-03
%   add bScaleAmp, if true, then scale amplitude of B s.t.
%   mean(B.amp(bMask)) = mean(A.amp(bMask))

sOptions = ValidateOptions(varargin{:});

% use function pointer to choose the desired AutoThreshold method
switch lower(sOptions.image_type),
    case 'pupil'
        AutoThreshold = @AutoThresholdPupil;
        
    case {'psf', 'source'}
        AutoThreshold = @AutoThresholdPSF;
        
    otherwise
        error(['unknown image type: ' sOptions.image_type]);
end % switch image_type

% if only one field is input, autothreshold and return
if ~exist('B','var') || isempty(B),
    [bMask, thresh] = AutoThreshold(A, sOptions);
    sResults = struct('thresh',thresh);
    return
end
   
% if two fields, first compare if they cover the same pixels
% create a mask if none given
if isempty(sOptions.Wgt_mask),
    
    [bMask_A, thresh] = AutoThreshold(A, sOptions);
    [bMask_B, thresh] = AutoThreshold(B, sOptions);
    
    bPassMask = isequal(bMask_A, bMask_B);

    bMask = bMask_A;
else
    
    % use the given mask, and no need to test
    [bMask, thresh] = AutoThreshold(sOptions.Wgt_mask);
    bPassMask = true;
    
end

% if allowed to adjust amplitude:
alpha = 1;
if sOptions.bScaleAmp,
    alpha = mean(abs(A(bMask)))./mean(abs(B(bMask)));
    B = alpha * B;
    if sOptions.debug,
        fprintf('Scale Amp, alpha = %e\n',alpha);
    end
end
    
% apply the metric in both the mask region...
% mf = max( |A-B| / min(|A|,|B|) )
mf_inside  = max( abs(A(bMask)-B(bMask)) ./ min(abs(A(bMask)),abs(B(bMask))) );

% and outside = ~bMask
% mf = max( |A-B| / mean(A(bMask)) )
mf_outside = max( abs(A(~bMask)-B(~bMask)) ./ mean(abs(A(bMask))) );

sResults = struct(...
    'bPassMask', bPassMask ...
    ,'bPassInside',  mf_inside  <= sOptions.eps_inside ...
    ,'bPassOutside', mf_outside <= sOptions.eps_outside ...
    ,'mf_inside', mf_inside ...
    ,'mf_outside', mf_outside ...
    ,'thresh', thresh ...
    ,'alpha', alpha ...
    );


end % main

%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bMask, thresh] = AutoThresholdPupil(Im, sOptions)
    % output bMask is boolean
    

    % use histogram to choose threshold
    [cnts, xbin] = hist(abs(Im(:)), sOptions.AutoThreshold_Nbins);
    [minima, xmin] = findpeaks(-cnts, xbin);

    function [amax, xmax] = findpeaks(y, x)
        
        bPeaks = [false (y(2:end-1) >= y(1:end-2) & y(2:end-1) > y(3:end)) false];
        
        amax = y(bPeaks);
        xmax = x(bPeaks);
        
    end % findpeaks

    % if there are no minima, minima is empty
    if isempty(minima),
        warning('could not find threshold, attempting AutoThresholdPSF');
        [bMask, thresh] = AutoThresholdPSF(Im, sOptions);
        return
    end

    % choose the first minimum
    thresh = xmin(1);

    bMask = abs(Im) > thresh;

    if sOptions.debug,
        figure, bar(xbin, cnts)
        hold on
        plot([1 1]'*xmin(:)', get(gca,'ylim')'*ones(1,length(xmin)),'-r')
        hold off
        grid on
    
        title(['Threshold = ' num2str(thresh)])
    end

end % AutoThreshold

function [bMask, thresh] = AutoThresholdPSF(Im, sOptions)
    % sOptions:
    %   logPSF
    %   PSF_thresh_nsig
    %   debug
    
    if sOptions.logPSF,
        absIm = abs(Im);
        %abmin = min(absIm(absIm(:)>0));
        absIm = absIm(absIm>0); % eliminate no signal
        [cnt, xbin] = hist(log10(absIm(:)), 256);
        xbin = 10.^xbin;
        
    else,
        [cnt, xbin] = hist(abs(Im(:)), 256);
        
    end
    
    % if debug, show amp and histogram
    if sOptions.debug,
        figure_mxn(1,2)
        subplot(1,2,1), imageschcit(real(log10(Im))), colorbartitle('log_{10} Im')
        subplot(1,2,2), plot(xbin, cnt,'-o'), grid
    end
    
    % assume most pixels are dark + noise
    [cntmax, imax] = max(cnt);
    
    % only pixels where Intensity is < twice the dark intensity
    Idk = Im(abs(Im) < 2*xbin(imax)); % xbin(imax+10)); 
    %[cntdk, xbindk] = hist(abs(Idk(:)), 256);
        
    % mean, std of dark pixels
    %     mu = sum(xbindk.*cntdk)./sum(cntdk);
    %     sig = sqrt( sum(xbindk.^2.*cntdk)./sum(cntdk) - mu.^2 );
    mu = mean(Idk(:));
    sig = std(Idk(:));
    
    thresh = mu + sOptions.PSF_thresh_nsig * sig;
    
    bMask = abs(Im) >= thresh;
    
    if sOptions.debug,
        fprintf('AutoThresholdPSF; thresh = %.3e\n',thresh);
        %figure, imageschcit(bMask), title('before imopen');
    end
    
    % morphological opening to eliminate outlier 'salt' noise
    bMask = imclose(bMask, strel('disk',3));
    bMask = imopen(bMask, strel('disk',3));

end % AutoThresholdPSF

function sOptions = ValidateOptions(varargin)
    
    % default options
    sOptions = struct(...
        'Wgt_mask', [] ... % option to supply a real valued array to create the mask
        ,'image_type', 'pupil' ... % 'pupil' or 'psf'
        ,'eps_inside', 1.0e-15 ...
        ,'eps_outside', 1.0e-12 ...
        ,'AutoThreshold_Nbins', 21 ...
        ,'bScaleAmp', false ...
        ,'logPSF', false ...
        ,'debug', false ...
        ,'PSF_thresh_nsig', 2 ...
        );
    
    % copy given option values over the default values
    if isempty(varargin), return, end

    % varargin can be a struct or list of keyword, value
    if isstruct(varargin{1})
        sTmp = varargin{1};
        fnames = fieldnames(sTmp);
        
        for ii = 1:length(fnames),
            sOptions.(fnames{ii}) = sTmp.(fnames{ii});
        end
        
    elseif iscell(varargin)
        fnames = fieldnames(sOptions);
        
        for ii = 1:length(fnames),
            fname = fnames{ii};
            sOptions.(fname) = CheckOption(fname, sOptions.(fname), varargin{:});            
        end

    end
    
    % validate option values here
    
        
end

