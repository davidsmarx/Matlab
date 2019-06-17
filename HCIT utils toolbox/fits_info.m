function header=fits_info(fitsfile)
%
% mfits: A simple fits package fo MATLAB.
%
% Author: Jason D. Fiege, University of Manitoba, 2010
% fiege@physics.umanitoba.ca
%
% The main component is a FITS writer fits_write to complement MATLAB's
% built-in fits reader.  fits_read is just a wrapper for fitsread.  fits_info
% returns the header information in a more convenient form than the
% built-in fitsinfo function.  Additional functionality is provides by
% getFitsValue and setFitsValue to get and set a single header field.
%
% -------------------------
% Just a wrapper for fitsinfo that returns the PrimaryData.Keywords
% part of the fits file as a MATLAB structure.  Only Keyword/value pairs
% are retained in the structure.  Comments are stripped.

% Read the fits file.  This returns an N x 2 cell array.
hCell=fitsinfo(fitsfile);
hCell=hCell.PrimaryData.Keywords(:,1:2);

% Convert to a structure
header=[];
for i=1:size(hCell, 1)
    keyword=hCell{i, 1};
    index=regexpi(keyword, '\W');
    keyword(index)=[];
    if ~isempty(keyword)
        header.(keyword)=hCell{i, 2};
    end
end

% Delete the 'END' field.
if isfield(header, 'END')
    header=rmfield(header, 'END');
end