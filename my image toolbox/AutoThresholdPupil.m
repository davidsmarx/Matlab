function [bMask, thresh] = AutoThresholdPupil(Im, varargin)
% [bMask, thresh] = AutoThresholdPupil(Im, sOptions)
% 
% find a threshold value, then apply the threshold to the image to
% create a binary mask
%
% the threshold is determined from the histogram of intensity values of Im
%
% Im = a 2-d image of real, positive values
% sOptions = struct( ...
%    'AutoThreshold_Nbins', (default = 21)
%    
% output bMask is boolean size of Im
    
% check options
sOptions = ValidateOptions(varargin{:});


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

function sOptions = ValidateOptions(varargin)
    
    % default options
    sOptions = struct(...
        'AutoThreshold_Nbins', 21 ...
        ,'debug', false ...
        );
    
    % check input for options
    if length(varargin) >= 1 && isstruct(varargin{1}),
        sOptin = varargin{1};
        fnames = fieldnames(sOptin);
        
        for ii = 1:length(fnames),
            sOptions.(fnames{ii}) = sTmp.(fnames{ii});
        end
    end % if input options
    
end % ValidateOptions


