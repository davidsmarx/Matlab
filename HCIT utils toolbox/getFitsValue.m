function value=getFitsValue(info, keyword)
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
if isfield(info,keyword)
    value=info.(keyword);
    if isnumeric(info.(keyword))
        value=double(value);
    end
else
    value=[];
end