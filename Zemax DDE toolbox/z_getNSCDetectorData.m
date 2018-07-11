function DetData = z_getNSCDetectorData(zchan, nsurf, nobject, pixel, datatype)
% DetData = z_getNSCDetectorData(zchan, nsurf, nobject, pixel, datatype)
%
% z_getNSCDetectorData(zchan) clears all detectors. Detectors must be
%       cleared before tracing and calling z_getNSCDetectorData() for retrieving
%       data
% 
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC, if nobject = 0, then all detectors are
%           cleared
% npixel = ix# is a positive integer
% greater than zero, then the data from the specified pixel is returned. Otherwise, the following data
% is returned depending upon the value of Pix#:
% 0: Sum of the total flux in position space, average flux/area in position space, or total flux in angle
% space for all pixels for that detector object, for Data = 0, 1, or 2, respectively.
% -1: Maximum flux or flux/area.
% -2: Minimum flux or flux/area.
% -3: Number of rays striking the detector.
% -4: Standard deviation (RMS from the mean) of all the non-zero pixel data.
% -5: The mean value of all the non-zero pixel data.
% -6, -7, -8: The x, y, or z coordinate of the position or angle Irradiance or Intensity centroid,
% respectively.
% -9,-10, -11, -12, -13: The RMS radius, x, y, z, or xy cross term distance or angle of all the pixel data
% with respect to the centroid. These are the second moments r^2, x^2, y^2, z^2, and xy, respectively.
%
% Datatype is 0 for flux, 1 for flux/area, 2 for flux/solid angle pixel, and 3 for normalized flux.

%     surf = ',1';
%     object = sprintf(',%d', sDetectorParms.DetectorObject);
%     pixel = ',0'; % return sum of all pixels for the detector
%     data = ',3'; % return power
% NSCDetectorData,surface,object,pixel,data

if nargin == 1,
    % clear all detectors
    cmdstr = 'NSCDetectorData,1,0,0,0';
    DetData = ddereq(zchan,cmdstr,[1 1]);

else
    
    cmdstr = ['NSCDetectorData,' num2str(nsurf)...
        ',' num2str(nobject)...
        ',' num2str(pixel)...
        ',' datatype];
    
    stmp = ddereq(zchan,cmdstr,[1 1]);
    
    DetData = str2double(stmp);
end
