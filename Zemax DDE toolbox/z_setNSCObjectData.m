function opval = z_setNSCObjectData(zchan, nsurf, nobject, sData)
% opval = z_setNSCObjectData(zchan, nsurf, nobject, sData)
% sData = z_setNSCObjectData
%    returns sData = struct with fields for all the possible data items
%
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC
% sData = struct of data items, see GetNSCObjectData, p. 826 for the table
%
% SetNSCObjectData, surf, object, code, data
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

if nargin == 0,
    opval = DefaultDataStruct;
    return
end

U = CConstants;

% for each non-empty field, look up the code and set the object value
if ~isstruct(sData),
    error('input sData is not a struct, use default struct');
end
opval = '';
fnames = fieldnames(sData);
for idata = 1:length(fnames),
    if isempty(sData.(fnames{idata})),
        continue;
    end
    
    % the data value
    data = sData.(fnames{idata});
    
    % look up the code and check valid data:
    emsg = '';
    strdata = ''; % the string passed with the DDE call
    switch fnames{idata},
        case 'ObjectTypeName',
            code = '0';
            if ~ischar(data), emsg = 'Type Name is a string'; end
            strdata = data;
        case 'Comment',
            code = '1';
            if ~ischar(data), emsg = 'Comment is a string'; end
            strdata = data;
        case 'Color',
            code = '2';
            strdata = num2str(data);
        case 'Reference', % object number. (integer)
            code = '5';
            strdata = num2str(data);
        case 'InsideOf', % [] ... object number. (integer)
            code = '6';
            strdata = num2str(data);
        case 'UserDefinedApertureFile', % [] ... 1 if object uses a user defined aperture file, 0 otherwise. (integer)
            code = '3';
            if data == 0 || data == 1,
                strdata = num2str(data);
            else,
                emsg = 'UserDefinedApertureFile must be 0 or 1';
            end
        case 'UserDefinedApertureFilename', % '' ... , if any. (string)
            code = '4';
            if ~ischar(data), emsg = 'UserDefinedApertureFilename must be string'; end
            strdata = data;
        case 'UsePixelInterpolation', % [] ... Sets the “Use Pixel Interpolation” checkbox. Use 1 for checked, 0 for unchecked.
            code = '29';
            strdata = num2str(data);
        case 'SourceRandomPolarization', % [] ... Sets the source object random polarization. Use 1 for checked, 0 for unchecked.
            code = '101';
            if data == 0 || data == 1,
                strdata = num2str(data);
            else
                emsg = 'SourceRandomPolarization must be 0 or 1';
            end
        case 'SourceReverseRays', % [] ... 102 Sets the source object reverse rays option. Use 1 for checked, 0 for unchecked.
            code = '102';
            if data == 0 || data == 1,
                strdata = num2str(data);
            else
                emsg = 'SourceReverseRays must be 0 or 1';
            end
        case 'SourceJonesX', % [] ... 103 Sets the source object Jones X value.
            code = '103';
            strdata = num2str(data);
        case 'SourceJonesY', % [] ... 104 Sets the source object Jones Y value.
            code = '104';
            strdata = num2str(data);
        case 'SourcePhaseX', % [] ... 105 Sets the source object Phase X value.
            code = '105';
            strdata = num2str(data/U.P); % in degrees
        case 'SourcePhaseY', % [] ... 106 Sets the source object Phase Y value.
            code = '106';
            strdata = num2str(data/U.P); % in degrees
        case 'SourceInitialPhase', % [] ...107 Sets the source object initial phase in degrees value.
            code = '107';
            strdata = num2str(data/U.P); % in degrees
        case 'SourceCoherenceLength', % [] ... 108 Sets the source object coherence length value.
            code = '108';
            strdata = num2str(data/U.MM); % in mm (lens units)
        case 'SourcePrePropagation', % [] ... 109 Sets the source object pre-propagation value.
            code = '109';
            strdata = num2str(data/U.MM); % in mm (lens units)
        case 'SourceSamplingMethod', % [] ... 110 Sets the source object sampling method; 0 for random, 1 for Sobol sampling.
            code = '110';
            if data == 0 || data == 1,
                strdata = num2str(data);
            else
                emsg = 'SourceSamplingMethod must be 0 or 1';
            end
        case 'SourceBulkScatter', % [] ... 111 Sets the source object bulk scatter method; 0 for many, 1 for once, 2 for never.
            code = '111';
            strdata = num2str(data);
        case 'BulkScatterMeanPathValue', % [] ... 202 Sets the Mean Path value.
            code = '202';
            strdata = num2str(data/U.MM); % mm
        case 'BulkScatterAngle', % [] ... 203 Sets the Angle value.
            code = '203';
            strdata = num2str(data/U.P); % in degrees
        case 'BulkScatterDLLParm', % [1-16, value] ... 211-226 Sets the DLL parameter 1-16, respectively.
            code = num2str(210 + data(1));
            strdata = num2str(data(2));
            
    end % switch fieldname
    
    % check for error
    if ~isempty(emsg),
        break;
    end
 
    % if strdata is empty, something went wrong
    if isempty(strdata) || isempty(code),
        warning(['property ' fnames{idata} ' has error']);
        continue;
    end
    
    % now send the command to ZEMAX
    % SetNSCObjectData, surf, object, code, data
     cmdstr = sprintf('SetNSCObjectData,%d,%d,%s,%s', ...
        nsurf,nobject,code,strdata);
    opval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

end % for each object property (field)

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
