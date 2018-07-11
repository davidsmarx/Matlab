function sData = z_getNSCObjectData(zchan, nsurf, nobject)
% sData = z_getNSCObjectData(zchan, nsurf, nobject)
% sData = z_getNSCObjectData
%    returns sData = struct with fields for all the possible data items
%
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC
% sData = struct of data items, see GetNSCObjectData, p. 826 for the table
%
% GetNSCObjectData, surf, object, code
% sData = struct(...
%   'ObjectTypeName', '' ... . (string)
%   ,'Comment', '' ... , which also defines the file name if the object is defined by a file. (string)
%   ,'Color', [] ... . (integer)
%   ,'Reference', [] ... object number. (integer)
%   ,'InsideOf', [] ... object number. (integer)
%   ... % The following codes set values on the Type tab of the Object Properties dialog.
%   ,'UserDefinedApertureFile', [] ... 1 if object uses a user defined aperture file, 0 otherwise. (integer)
%   ,'UserDefinedApertureFilename', '' ... , if any. (string)
%   ,'UsePixelInterpolation', [] ... Sets the “Use Pixel Interpolation” checkbox. Use 1 for checked, 0 for unchecked.
%   ... % The following codes set values on the Sources tab of the Object Properties dialog.
%   ,'SourceRandomPolarization', [] ... Sets the source object random polarization. Use 1 for checked, 0 for unchecked.
%   ,'SourceReverseRays', [] ... 102 Sets the source object reverse rays option. Use 1 for checked, 0 for unchecked.
%   ,'SourceJonesX', [] ... 103 Sets the source object Jones X value.
%   ,'SourceJonesY', [] ... 104 Sets the source object Jones Y value.
%   ,'SourcePhaseX', [] ... 105 Sets the source object Phase X value.
%   ,'SourcePhaseY', [] ... 106 Sets the source object Phase Y value.
%   ,'SourceInitialPhase', [] ...107 Sets the source object initial phase in degrees value.
%   ,'SourceCoherenceLength', [] ... 108 Sets the source object coherence length value.
%   ,'SourcePrePropagation', [] ... 109 Sets the source object pre-propagation value.
%   ,'SourceSamplingMethod', [] ... 110 Sets the source object sampling method; 0 for random, 1 for Sobol sampling.
%   ,'SourceBulkScatter', [] ... 111 Sets the source object bulk scatter method; 0 for many, 1 for once, 2 for never.
%   ... % The following codes set values on the Bulk Scatter tab of the Object Properties dialog.
%   ,'BulkScatterMeanPathValue', [] ... 202 Sets the Mean Path value.
%   ,'BulkScatterAngle', [] ... 203 Sets the Angle value.
%   ,'BulkScatterDLLParm', [[],[]] ... 211-226 Sets the DLL parameter 1-16, respectively.
%   );


% Code Data returned by GetNSCObjectData
%   0 Object type name. (string)
%   1 Comment, which also defines the file name if the object is defined by a file. (string)
%   2 Color. (integer)
%   5 Reference object number. (integer)
%   6 Inside of object number. (integer)
%     The following codes set values on the Type tab of the Object Properties dialog.
%   3 1 if object uses a user defined aperture file, 0 otherwise. (integer)
%   4 User defined aperture file name, if any. (string)
%  29 Sets the “Use Pixel Interpolation” checkbox. Use 1 for checked, 0 for unchecked.
%     The following codes set values on the Sources tab of the Object Properties dialog.
% 101 Sets the source object random polarization. Use 1 for checked, 0 for unchecked.
% 102 Sets the source object reverse rays option. Use 1 for checked, 0 for unchecked.
% 103 Sets the source object Jones X value.
% 104 Sets the source object Jones Y value.
% 105 Sets the source object Phase X value.
% 106 Sets the source object Phase Y value.
% 107 Sets the source object initial phase in degrees value.
% 108 Sets the source object coherence length value.
% 109 Sets the source object pre-propagation value.
% 110 Sets the source object sampling method; 0 for random, 1 for Sobol sampling.
% 111 Sets the source object bulk scatter method; 0 for many, 1 for once, 2 for never.
%     The following codes set values on the Bulk Scatter tab of the Object Properties dialog.
% 202 Sets the Mean Path value.
% 203 Sets the Angle value.
% 211-226 Sets the DLL parameter 1-16, respectively.

U = CConstants;

sData = DefaultDataStruct;

if nargin == 0,
    return
end    
    
sData.ObjectTypeName = GetData('0','str');
sData.Comment = GetData('1','str');
sData.Color = GetData('2',1);
sData.Reference = GetData('5',1);
sData.InsideOf = GetData('6',1);
sData.UserDefinedApertureFile = GetData('3',1); % 1 = on, 0 = off
sData.UserDefinedApertureFilename = GetData('4','str');
sData.UsePixelInterpolation = GetData('29',1);
sData.SourceRandomPolarization = GetData('101',1);
sData.SourceReverseRays = GetData('102',1);
sData.SourceJonesX = GetData('103',1.0); % normalized field value
sData.SourceJonesY = GetData('104',1.0);
sData.SourcePhaseX = GetData('105', U.P); % phase in degrees
sData.SourcePhaseY = GetData('106', U.P); % phase in degrees
sData.SourceInitialPhase = GetData('107', U.P); % phase in degrees
sData.SourceCoherenceLength = GetData('108', U.MM); % lens units
sData.SourcePrePropagation = GetData('109', U.MM); % lens units
sData.SourceSamplingMethod = GetData('110',1);
sData.SourceBulkScatter = GetData('111',1);
sData.BulkScatterMeanPathValue = GetData('202',U.MM);
sData.BulkScatterAngle = GetData('203',U.P);
sData.BulkScatterDLLParm = []; % not implemented



    function value = GetData(code,sUnits)

        % now send the command to ZEMAX
        % GetNSCObjectData, surf, object, code
        cmdstr = sprintf('GetNSCObjectData,%d,%d,%s', ...
            nsurf,nobject,code);
        opvalstr = ddereq(zchan,cmdstr,[1 1]);
        
        if ischar(sUnits),
            value = opvalstr;
        else
            value = str2double(opvalstr)*sUnits;
        end
            
    end % function GetData

end % main


function sData = DefaultDataStruct


sData = struct(...
  'ObjectTypeName', '' ... . (string)
  ,'Comment', '' ... , which also defines the file name if the object is defined by a file. (string)
  ,'Color', [] ... . (integer)
  ,'Reference', [] ... object number. (integer)
  ,'InsideOf', [] ... object number. (integer)
  ... % The following codes set values on the Type tab of the Object Properties dialog.
  ,'UserDefinedApertureFile', [] ... 1 if object uses a user defined aperture file, 0 otherwise. (integer)
  ,'UserDefinedApertureFilename', '' ... , if any. (string)
  ,'UsePixelInterpolation', [] ... Sets the “Use Pixel Interpolation” checkbox. Use 1 for checked, 0 for unchecked.
  ... % The following codes set values on the Sources tab of the Object Properties dialog.
  ,'SourceRandomPolarization', [] ... Sets the source object random polarization. Use 1 for checked, 0 for unchecked.
  ,'SourceReverseRays', [] ... 102 Sets the source object reverse rays option. Use 1 for checked, 0 for unchecked.
  ,'SourceJonesX', [] ... 103 Sets the source object Jones X value.
  ,'SourceJonesY', [] ... 104 Sets the source object Jones Y value.
  ,'SourcePhaseX', [] ... 105 Sets the source object Phase X value.
  ,'SourcePhaseY', [] ... 106 Sets the source object Phase Y value.
  ,'SourceInitialPhase', [] ...107 Sets the source object initial phase in degrees value.
  ,'SourceCoherenceLength', [] ... 108 Sets the source object coherence length value.
  ,'SourcePrePropagation', [] ... 109 Sets the source object pre-propagation value.
  ,'SourceSamplingMethod', [] ... 110 Sets the source object sampling method; 0 for random, 1 for Sobol sampling.
  ,'SourceBulkScatter', [] ... 111 Sets the source object bulk scatter method; 0 for many, 1 for once, 2 for never.
  ... % The following codes set values on the Bulk Scatter tab of the Object Properties dialog.
  ,'BulkScatterMeanPathValue', [] ... 202 Sets the Mean Path value.
  ,'BulkScatterAngle', [] ... 203 Sets the Angle value.
  ,'BulkScatterDLLParm', [[],[]] ... 211-226 Sets the DLL parameter 1-16, respectively.
  );

end % DefaultDataStruct
