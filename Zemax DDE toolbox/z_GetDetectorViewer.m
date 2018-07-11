function [Data, sParms] = z_GetDetectorViewer(varargin)
% [Data, sParms] = z_GetDetectorViewer(zchan, textfilename)
% 
% [Data, sParms] = z_GetDetectorViewer(zchan)
% [Data, sParms] = z_GetDetectorViewer(textfilename)
% [Data, sParms] = z_GetDetectorViewer(zchan, textfilename)
% [Data, sParms] = z_GetDetectorViewer(zchan, textfilename, settings.cfg)
% [Data, sParms] = z_GetDetectorViewer
%
% for NSC systems
% calls ZEMAX Analysis -> Ray Tracing -> Detector Dataer
% the settings are the default settings for the lens
%
% [...] = z_GetDetectorViewer(zchan)
%    ZEMAX calculates the Datamap on current lens, results are stored in
%    [pwd '\detectorviewer.txt']
%
% [...] = z_GetDetectorViewer(textfilename)
%    read results from textfilename
%
% [...] = z_GetDetectorViewer(zchan, textfilename)
%    ZEMAX calculates the wavefrontmap on current lens and stores results in
%    textfilename
% 
% [...] = z_GetDetectorViewer
%    interactively select a text file with wavefrontmap results to read
%
% [...] = z_GetDetectorViewer(zchan, textfilename, settings.cfg)
%    use the configuration in settings.cfg to make the calculations. If
%    settings.cfg is not specified, default settings are used, and then the
%    default settings are written to 
%
% Data = pixelated detector view
% sParms fields:
%    file
%    date
%    title
%    configuration
%    description
%    wavelength
%    field
%    PeakToValley, RMS
%    surface
%    exitpupildiam
%    GridSize [x,y] (pixels)
%    GridCenter [x0, y0] (pixels)

% validate inputs
[zchan, textfilename, settingsfilename, settingsflag] =...
    ValidateInputs(varargin{:});

% if zchan is empty, then just want to read an existing text file
if ~isempty(zchan),
    % Geometric Image Analysis
    typecode = 'Dvr';
    
    % use temp file for textfile to keep cmdstr length < 255
    tempfilename = 'c:\temp\dvrtextfile.txt';
    
    % make the call
    cmdstr =...
        ['GetTextFile, "' tempfilename '", ' typecode ', "' settingsfilename '", ' settingsflag];
    if length(cmdstr) >= 255, error(['cmdstr length = ' num2str(length(cmdstr))]); end
    retval = ddereq(zchan,cmdstr,[1 1]);
    
    % wait until zemax is ready
    while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
    
    % now copy temp file to specified file
    copyfile(tempfilename, textfilename);

end

% get the image data from the textfile
[Data, sParms] = ReadDatamap(textfilename);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Data, sParms] = ReadDatamap(textfile)

U = CConstants;

% sParms fields:
%    file
%    date
%    configuration
%    description
%    wavelength
%    wavelength
%    field
%    PeakToValley, RMS
%    surface
%    exitpupildiam
%    GridSize [x,y] (pixels)
%    GridCenter [x0, y0] (pixels)

% initialize outputs
Data = [];
sParms = BlankParmsStruct;

fid = fopen(textfile,'rt');
if fid == -1, error(['error opening text file ' textfile]); end

try,
    
    ltmp = fgetl(fid);

while ~feof(fid)
    if ~isempty(ltmp)
        wordparse = textscan(ltmp,'%s');

        switch wordparse{1}{1}, % first word
            case 'File',
                atmp = textscan(ltmp,'File : %s','whitespace','\t');
                sParms.Lensfile = atmp{1}{1};
            case 'Title:',
                atmp = textscan(ltmp,'Title: %s','whitespace','\t');
                if ~isempty(atmp{1}), sParms.Title = atmp{1}{1}; end
            case 'Date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                sParms.Date = atmp{1}{1};
            case 'Configuration'
                atmp = textscan(ltmp,'Configuration %d of %d');
                sParms.Configuration = atmp{1};
            case 'Detector'
                switch wordparse{1}{2}
                    case {'X', 'Y', 'Z'},
                        atmp = textscan(wordparse{1}{end},'%f');
                        sParms.(wordparse{1}{2}) = atmp{1}*U.MM;

                    case 'Tilt'
                        atmp = textscan(wordparse{1}{end},'%f');
                        sParms.([wordparse{1}{2} '_' wordparse{1}{3}]) = atmp{1}*U.P;
                        
                    case 'Viewer'
                        % do nothing
                        
                    otherwise
                        atmp = textscan(wordparse{1}{2},'%d,');
                        sParms.DetectorObject = atmp{1};

                end % switch Detector                   

            case 'Size'
                atmp = textscan(ltmp,'Size %f W X %f H %s Pixels %f W X %f H Total Hits = %f','delimiter',',');
                sParms.DetectorSize = [atmp{1} atmp{2}]*U.MM;
                sParms.DetectorPixels = [atmp{4} atmp{5}];
                sParms.TotalHits = atmp{end};
                
                
            case 'Peak' % Irradiance : 1.6459E+000 Watts/cm^2
                atmp = textscan(ltmp,'Peak Irradiance : %f');
                sParms.PeakIrradiance = atmp{1};
                
            case 'Total' % Power     : 5.1619E-002 Watts
                atmp = textscan(ltmp,'Total Power     : %f');
                sParms.TotalPower = atmp{1};

            case 'Smoothing'
                atmp = textscan(ltmp,'Smoothing : %s');
                sParms.Smoothing = atmp{1}{1};
            case 'Data'
                atmp = textscan(ltmp,'Data Type :  %s','whitespace','\t');
                sParms.DataType = atmp{1}{1};
            
            case 'Position'
                atmp = textscan(ltmp,'Position Units : %s');
                
            case 'Units'
                atmp = textscan(ltmp,'Units : %s');
                dataunits = atmp{1}{1};
                
                % data follows
                [Data, sParms] = ReadData(fid, sParms);
                
                
            otherwise,

        end % switch

    end % while ~isempty
    
    % get next line
    ltmp = fgetl(fid);

end % while ~feof

switch dataunits
    case 'degrees'
        Data = Data * U.P;
    case 'Watts/cm^2';
        % do nothing
        
    otherwise
        error(['unknown units: ' dataunits]);
end

catch
    fprintf('encountered error, closing file\n');
    
end % try

fclose(fid);

end % main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Data, sParms] = ReadData(fid, sParms)

% skip empty lines
ltmp = SkipEmptyLines(fid);

% check if profile or 2-d image map
Atmp = textscan(ltmp,'%s');

switch Atmp{1}{1}
    
    case 'Row'
        ltmp = SkipEmptyLines(fid);
        % ltmp is now = 'X Coordinate    Value', (x,v) data follows
        Data = fscanf(fid,'%f %f',[2 sParms.DetectorPixels(2)])';
        
    case 'Column'
        ltmp = SkipEmptyLines(fid);
        % ltmp is now = 'Y Coordinate    Value', (x,v) data follows
        Data = fscanf(fid,'%f %f',[2 sParms.DetectorPixels(1)])';

    case '1' % 2-d array type, header line of column #'s

        dtmp = fscanf(fid,'%f',...
            [sParms.DetectorPixels(1)+1 sParms.DetectorPixels(2)]);
        Data = dtmp(2:end,:)';

        % flip Data array in y-direction because the first data point in
        % the text file is
        % the upper left corner, which is (-x, +y), but in the Matlab array
        % the first point is (-x, -y) in the grid
        Data = flipud(Data);
        
    otherwise,
        Data = [];
        
end % switch data type


end % ReadData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ltmp = SkipEmptyLines(fid)

    ltmp = [];
    while isempty(ltmp) && ~feof(fid),
        ltmp = fgetl(fid);
    end

end % SkipEmptyLines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sParms = BlankParmsStruct

    sParms = struct('Lensfile', [] ...
        , 'Title', [] ...
        , 'Date', [] ...
        , 'DetectorObject', [] ...
        , 'DetectorSize', [] ...
        , 'DetectorPixels', [] ...
        , 'TotalHits', [] ...
        , 'PeakIrradiance', [] ...
        , 'TotalPower' , [] ...
        , 'Smoothing', [] ...
        , 'DataType', [] ...
        , 'X', [] ...
        , 'Y', [] ...
        , 'Z', [] ...
        , 'Tilt_X', [] ...
        , 'Tilt_Y', [] ...
        , 'Tilt_Z', [] ...
        );
    
end % BlankParmsStruct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin)
% [Field, sParms] = z_GetDetectorViewer(zchan)
% [Field, sParms] = z_GetDetectorViewer(textfilename)
% [Field, sParms] = z_GetDetectorViewer(zchan, textfilename)
% [Field, sParms] = z_GetDetectorViewer(zchan, textfilename, settings.cfg)

% defaults:
zchan = [];
settingsfilename = [pwd '\DvrSettings.cfg'];
textfilename = [pwd '\dvranalysis.txt'];
settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename

switch(nargin),
    case 0,
        [filename, pathname, fi] = uigetfile({'*.txt'},'Select wavefrontmap Analysis');
        if ~isequal(filename,0),
            textfilename = [pathname filename];
        end

        % use defaults for everything else

    case 1,
        if isstr(varargin{1}),
            textfilename = varargin{1};
        else
            zchan = varargin{1};
        end

        % use defaults for settings
        
    case 2,
        zchan = varargin{1};
        textfilename = varargin{2};
        
    case 3,
        [zchan, textfilename, settingsfilename] = deal(varargin{:});
        settingsflag = '1'; % 1 => use settings in settingsfilename

    otherwise,
        error('usage: [Field, sParms] = z_GetDetectorViewer(zchan, textfilename)');
end

end % ValidateInputs
