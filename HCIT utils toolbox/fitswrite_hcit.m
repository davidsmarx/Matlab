function fitswrite_hcit(out_fn, img, kwds, varargin)
% fitswrite_hcit(fn, img, kwds, 'image', {img, kwds}, ...)
% 
% my fits write wrapper to include headers and multiple hdus
%
% kwds is n x 3 cell array, just like finfo.PrimaryData.Keywords
% options:
%   'overwrite', false (default), or true

import matlab.io.*

% options
bOverwrite = CheckOption('overwrite', false, varargin{:});

if exist(out_fn, 'file') && bOverwrite,
    delete(out_fn);
elseif exist(out_fn, 'file') && ~bOverwrite,
    error(['file ' out_fn ' exists and overwrite is false']);    
end

% create output fits file
fptr = fits.createFile(out_fn);

try, % if error, close fptr
    
    % primary hdu
    fits.createImg(fptr,'double_img',size(img));
    fits.writeImg(fptr, img);
    [kwds, Nkeys] = CheckKwds(kwds);
    for ikey = 1:Nkeys,
        fits.writeKey(fptr, kwds{ikey,:});
    end
    
    % if additional image hdus
    ii_images = find(strcmp(varargin, 'image'));
    for ii_imagehdu = 1:length(ii_images),
        
        img = varargin{ii_images(ii_imagehdu)+1}{1};
        kwds = varargin{ii_images(ii_imagehdu)+1}{2};
        
        fits.createImg(fptr,'double_img',size(img))
        fits.writeImg(fptr, img);
        Nkeys = CheckKwds(kwds);
        for ikey = 1:Nkeys,
            fits.writeKey(fptr, kwds{ikey,:});
        end
        
    end % for each image hdu

catch ME
    warning('there was an error:');
    disp(ME);
end % try

fits.closeFile(fptr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [kwds, Nkeys] = CheckKwds(kwds)
        % check kwds is valid Nkeys x 3 cell array
        % remove rows with reserved words
        % remove rows with empty values
        
        if ~iscell(kwds), error(['kwds is not a cell array']); end
        
        if isempty(kwds), Nkeys = 0; return, end
        
        % list of reserved words
        list_reserved_keys = {
            'END'
            'SIMPLE'
            'BITPIX'
            };
        for ik = 1:length(list_reserved_keys)
            irow = strcmp(kwds(:,1), list_reserved_keys{ik});
            if ~isempty(irow)
                kwds(irow,:) = [];
            end
        end

        % reserved word NAXIS is special
        irow = strncmp(kwds(:,1), 'NAXIS', length('NAXIS'));
        if ~isempty(irow)
            kwds(irow,:) = [];
        end
        
        % rows with empty values
        irows = find(cellfun(@isempty,kwds(:,2)));
        irows = flipud(irows); % go in reverse order so location of empty rows doesn't change as rows are deleted
        for ik = 1:length(irows)
            kwds(irows(ik),:) = [];
        end
        
        % rows with empty comments, replace comment with a string
        irows = find(cellfun(@isempty,kwds(:,3)));
        for ik = 1:length(irows)
            kwds(irows(id),3) = ' ';
        end
        
        [Nkeys, n3] = size(kwds);
        if ~isequal(n3, 3),
            fits.closeFile(fptr);
            error(['kwds does not have 3 columns']);
        end
        
    end % CheckKwds

end % main