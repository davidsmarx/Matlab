function aout = grayto8(ain)
%GRAYTO8 Scale and convert grayscale image to uint8.
%   B = GRAYTO8(A) converts the double array A to uint8 by
%   scaling A by 255 and then rounding.  NaN's in A are converted
%   to 0.  Values in A greater than 1.0 are converted to 255;
%   values less than 0.0 are converted to 0.
%
%   B = GRAYTO8(A) converts the uint16 array A by scaling the
%   elements of A by 1/257, rounding, and then casting to uint8.
%
%   Copyright 1993-1998 The MathWorks, Inc. All Rights Reserved.
%   $Revision: 1.2 $  $Date: 1998/03/31 22:53:28 $

ain(isnan(ain)) = 0;

aout = uint8(round(255.*ain));

return
