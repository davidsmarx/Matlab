function [Field, headerinfo] = z_GetExtendedCohImage(varargin)
% [Img, parmsstruct] = z_ExtendedCohImage(zchan, options)
% [Img, parmsstruct] = z_ExtendedCohImage(options)
% 
% options:
%     'textfilename', textfilename (default = [pwd '\ExtCohImg.txt'])
%     'settingsfilename', settingsfilename ( .cfg ), (default = use lens
%                                                     current settings)
%
% calls ZEMAX Huygens PSF analysis
% the settings are the default settings for the lens
%
%
% Field = matrix of field values calculated from the
% headerinfo:

% validate inputs
[zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin{:});

if ~isempty(zchan),
    try
        % call ZEMAX to do the analysis
        
        % zemax three-letter code for this analysis
        % Xdi EXTENDED_DIFFRACTION_IMAGE_ANALYSIS
        typecode = 'Xdi';
        
        % make the call
        cmdstr =...
            ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', ' settingsflag];
        retval = ddereq(zchan,cmdstr,[1 1]);
        
        % wait until zemax is ready
        while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
    catch
        disp('zemax dde error!');
        keyboard;
    end % try catch
end

% get the image data from the textfile
[Field, headerinfo] = ReadAnalysisTxt(textfilename);
headerinfo.textfilename = textfilename;
headerinfo.settingsfilename = settingsfilename;
headerinfo.settingsflag = settingsflag;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Img, headerinfo] = ReadAnalysisTxt(textfilename)

global NM UM MM P;

fid = fopen(textfilename,'rt');
if fid == -1, error(['error opening text file ' textfilename]); end

% read header info
% headerinfo:

% Listing of Extended Diffraction Image Analysis Data
% 
% File : Y:\WFIRST\PupilImageQuality_20181128\WFIRST_20200105_Zemax\WFIRST-TCA-20200105_CGI-DI-20191025_TRANSLATION_W_OBSCURATION_4_PupilObject_to_SPAM_20200414.ZMX
% Title: WFIRST-TCA-20200105_CGI-DI-20191025.LEN
% Date : 4/17/2020
% 
% 
% 0.7300 ?m. Position: 0.00, 0.00 mm
% Source : 18.0000 Millimeters wide, 1002 by 1002 pixels
% Display: is 19.0000 Millimeters wide, 1024 by 1024 pixels
% Diffraction Limited Response
% Fractional transmission of image surface aperture: 1.000000

% 
%    Lensfile
%    Title
%    Date
%
%    Wavelength = [start_wave end_wave]
%    Field Position
%    Source dimension and sampling
%    Display dimension and sampling
%    "Diffraction Limited Response"
%    Fractional transmission of image surface aperture: 1.000000

% initialize return vals to empty, if something goes wrong with reading the
% text file, we return empty
Img = [];
headerinfo = struct;

% read each line, assumes 'Values...' is the last header line before data
% not all header lines are always included, depending on settings
ltmp = fgetl(fid); % 
ltmp = fgetl(fid);

while ~feof(fid),
    if ~isempty(ltmp),
        wordparse = textscan(ltmp,'%s');
    
        switch lower(wordparse{1}{1}), % first word
            case 'file',
                atmp = textscan(ltmp,'File : %s','whitespace','\t');
                headerinfo.Lensfile = atmp{1}{1};
            case 'title:',
                atmp = textscan(ltmp,'%s %s','delimiter',':');
                headerinfo.Title = atmp{2}{1};
            case 'date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                headerinfo.Date = atmp{1}{1};

                % a couple of blank lines, then wavelength and field point
                ltmp = fgetl(fid);
                ltmp = fgetl(fid); ltmp = fgetl(fid);
                ctmp = textscan(ltmp,'%f %s Position: %f, %f'); % % 0.7300 ?m. Position: 0.00, 0.00 mm
                headerinfo.wave = ctmp{1}*UM;
                headerinfo.fieldpoint = [ctmp{3:4}]*MM;
                headerinfo.fieldpointUnits = 'm';

            %case '', % wavelength, come back to this
            
            case 'source',
                % Source : 18.0000 Millimeters wide, 1002 by 1002 pixels
                ctmp = textscan(ltmp,'Source : %f Millimeters wide, %d by %d');
                headerinfo.SourceWidth = ctmp{1}*MM;
                headerinfo.SourceNxNy = [ctmp{2:3}];
                
            case 'display:'
                % Display: is 19.0000 Millimeters wide, 1024 by 1024 pixels
                ctmp = textscan(ltmp,'Display: is %f Millimeters wide, %d by %d');
                headerinfo.DisplayWidth = ctmp{1}*MM;
                headerinfo.DisplayNxNy = [ctmp{2:3}];
                
            case 'diffraction'
                headerinfo.DiffractionComment = ltmp;
                
            case 'fractional'
                headerinfo.Fractional = textscan(ltmp,'%s %f','delimiter',':');
                
                % a couple of blank lines, then data
                ltmp = fgetl(fid); % blank line before data
                ltmp = fgetl(fid); % blank line before data
                
                % read the data
                Img = fscanf(fid,'%e',headerinfo.DisplayNxNy)';
                
                %             otherwise,
                %                 % if line begins with a number, assume it is in the format:
                %                 % wavelength um at fieldpoint_x, fieldpoint_y units.
                %                 if ~isnan(str2double(wordparse{1}{1})),
                %                     headerinfo.Wavelegnth = str2double(wordparse{1}{1})*UM;
                %                     headerinfo.fieldpoint = str2double(wordparse{1}(4:5));
                %                     headerinfo.fieldpointUnits = wordparse{1}{end};
                %
                %                     if strfind(headerinfo.fieldpointUnits, 'mm')
                %                         headerinfo.fieldpoint = headerinfo.fieldpoint*MM;
                %                     elseif strfind(headerinfo.fieldpointUnits, 'deg')
                %                         headerinfo.fieldpoint = headerinfo.fieldpoint*P;
                %                     end
                %                 end
                
        end % switch
    end % if ~isempty line

    % read next line
    ltmp = fgetl(fid);

end % while ~feof

fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan, options)
% [Field, parmsstruct] = z_GetPOPAnalysis(options)

zchan = [];
if nargin >= 1 && isa(varargin{1},'double'),
    zchan = varargin{1};
end

textfilename = CheckOption('textfilename', [pwd '\PSFHuygens.txt'], varargin{:});
settingsfilename = CheckOption('settingsfilename', '', varargin{:});

if isempty(settingsfilename),
    settingsfilename = [pwd '\Settings.cfg'];
    settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename, see p. 634
else
    settingsflag = '1'; % 1 => use settings in settingsfilename
end

% switch(nargin),
%     case 0,
%         [filename, pathname, fi] = uigetfile({'*.txt'},'Select PSFHuygens Analysis');
%         if isequal(filename,0),
%             textfilename = [];
%         else,
%             textfilename = [pathname filename];
%         end
%         zchan = [];
% end

end % ValidateInputs
