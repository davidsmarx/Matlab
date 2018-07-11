function emsg = z_ModifySettings( zchan, cfgfilename, typecode, value )
% emsg = z_ModifySettings( ZEMAXDDE, cfgfilename, type, value )
%   
% type = string mnemonic, see manual
% some types from the manual:
%
% %%%% Detector Viewer: %%%%%
% DVW_SURFACE: The surface number. Use 1 for Non-Sequential mode.
% DVW_DETECTOR: The detector number.
%
% DVW_SHOW: The “show as” setting. The meaning depends upon the type of
% window displayed:
% For Graphics Windows: Use 0 for grey scale, 1 for inverted grey scale, 2 for false
% color, 3 for inverted false color, 4 for cross section row, and 5 for cross section
% column.
% For Text Windows: Use 0 for full pixel data, 1 for cross section row, and 2 for
% cross section column.
%
% DVW_ROWCOL: The row or column number for cross section plots.
%
% DVW_ZPLANE: The Z-Plane number for detector volumes.
%
% DVW_SCALE: The scale mode. Use 0 for linear, 1 for Log -5, 2 for Log -10, and
% 3 for Log - 15.
%
% DVW_SMOOTHING: The integer smoothing value.
%
% DVW_DATA: Use 0 for incoherent irradiance, 1 for coherent irradiance, 2 for
% coherent phase, 3 for radiant intensity, 4 for radiance (position space), and 5 for
% radiance (angle space).
%
% DVW_ZRD: The ray data base name, or null for none.
% DVW_FILTER: The filter string.
% DVW_MAXPLOT: The maximum plot scale.
% DVW_MINPLOT: The minimum plot scale.
% DVW_OUTPUTFILE: The output file name.
%
% %%%%% POP Analysis %%%%%
%
% POP_BEAMTYPE
%   0 => Gaussian Waist
%   4 => ZBF file
%
% For Beam Type = Gaussian Waist
% POP_PARAMn, where
%   n = 1 => Waist X
%   n = 2 => Waist Y
%   n = 3 => Decenter X
%   n = 4 => Decenter Y
%   n = 5 => Aperture X
%   n = 6 => Aperture Y
%   n = 7 => Order X
%   n = 8 => Order Y
%
% For Beam type = File
% POP_SOURCEFILE: name of ZBF file


emsg = '';

% ModifySettings, filename, type, value

% value should be a string, if not convert
if ~ischar(value),
    vtmp = value;
    value = num2str(vtmp, '%.5e');
end

cmdstr =...
    ['ModifySettings, "' cfgfilename '", ' typecode ', "' value '"'];
if length(cmdstr) >= 255, error(['cmdstr length = ' num2str(length(cmdstr))]); end

retval = ddereq(zchan,cmdstr,[1 1]);
if ~ischar(retval),
    disp(retval);
    retval = num2str(retval);
    %keyboard;
end

switch str2num(retval)
    case 0,
        emsg = ''; % no error
    case -1,
        emsg = 'Invalid File';
    case -2,
        emsg = 'Incorrect Version Number';
    case -3,
        emsg = 'File Access Conflict';
    otherwise,
        emsg = 'Unknown Error';
end


end

