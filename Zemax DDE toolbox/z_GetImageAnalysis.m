function [efficiency, img, field, pupil, wave] = z_GetImageAnalysis(varargin)
% Sout = z_GetImageAnalysis(zchan)
% Sout = z_GetImageAnalysis(textfilename)
% Sout = z_GetImageAnalysis(zchan, textfilename)
% Sout = z_GetImageAnalysis(zchan, textfilename, settings.cfg)
% [efficiency, img, field, pupil, wave] = z_GetImageAnalysis(zchan, textfilename)
% 
% calls ZEMAX geometric image analysis
% the settings are the default settings for the lens
%
% Results are stored in ASCII form in textfilename, the
% default is [pwd '\imageanalysis.txt']
% 
% efficiency = fraction of rays to reach the image plane
% img = [x,y] location of rays in the image plane
% field = [x,y] location of rays at source
% pupil = [x,y] location of rays at the pupil
% wave = wavelength# for each ray


[zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin{:});

if ~isempty(zchan),
    % Geometric Image Analysis
    typecode = 'Ima';
    
    % make the call
    cmdstr =...
        ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', ' settingsflag];
    retval = ddereq(zchan,cmdstr,[1 1]);
    
    % wait until zemax is ready
    while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
end

% get the image data from the textfile
sParms = ReadImageanalysisTxt(textfilename);

% add x,y vectors for convenience
pixelsize = sParms.ImageWidth ./ sParms.NumPixels;
firstpixel = -sParms.ImageWidth/2 + pixelsize/2;
lastpixel = -firstpixel;
sParms.x = linspace(firstpixel(1),lastpixel(1),sParms.NumPixels(1))';
sParms.y = linspace(firstpixel(2),lastpixel(2),sParms.NumPixels(2))';


% for backward compatibility efficiency, img, field, pupil, wave
if nargout > 1,
    field = sParms.field;
    pupil = sParms.field;
    wave  = sParms.wave;
    img   = sParms.img;
    efficiency = sParms.efficiency;
    
else,
    efficiency = sParms;

end

end % main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sParms = ReadImageanalysisTxt(textfile)

U = CConstants;

sParms = struct(...
    'field',[]...
    ,'pupil',[]...
    ,'wave',[]...
    ,'img',[]...
    ,'efficiency',[]...
    ,'Datatype', []...
    ,'ImageWidth', []...
    ,'NumPixels', []...
    ,'Units', []...
    ,'Lensfile',[]...
    ,'Title',[]...
    ,'Date',[]...
    ,'ObjectField',[]...
    );

fid = fopen(textfile,'rt');
if fid == -1, error(['error opening text file ' textfile]); end

try
    
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

            case 'Object',
                atmp = textscan(ltmp,'Object field: %f,%f');
                sParms.ObjectField = [atmp{1:2}];
                
            case 'Ray', % spot diagram type, ray listing follows
                sParms.Datatype = 'ray table';
                % data follows
                % read in all of the ray data
                C = textscan(fid,'%d %f %f %f %f %f %f %f %f','headerlines',12);
                [rayno, xfield, yfield, px, py, wave, wgt, xi, yi] = deal(C{:});
                sParms.field = [xfield yfield]*U.MM;
                sParms.pupil = [px py]*U.MM;
                sParms.img   = [xi yi]*U.MM;           
                sParms.wave  = wave;
                
            case 'Image' % Width                    : 20 Millimeters
                atmp = textscan(ltmp, 'Image Width : %f');
                sParms.ImageWidth = atmp{1}*U.MM;

            case 'Number', % of pixels
                atmp = textscan(ltmp,'Number of pixels : %f x %f');
                sParms.NumPixels = [atmp{1:2}];
                
            case 'Percent', % Percent Efficiency
                atmp = textscan(ltmp,'Percent Efficiency : %f');
                sParms.efficiency = atmp{1}/100;
                
            case 'Units', % histogram type, intensity per pixel follows
                sParms.Datatype = 'histogram';
                sParms.Units = textscan(ltmp,'Units : %s'); % Watts/millimeter squared
                                
                sParms.img = fscanf(fid,'%f', ...
                     [sParms.NumPixels(1) sParms.NumPixels(2)])'...
                     ./(U.MM*U.MM); % to W/m^2

            otherwise,

        end % switch

    end % while ~isempty
    
    % get next line
    ltmp = fgetl(fid);

end % while ~feof


% skip to last line to get efficiency value
if fseek(fid,-52,'eof') == -1, error(ferror(fid)); end
efficiency = fscanf(fid,'%*s %f');
efficiency = 0.01*efficiency; % convert percent to fraction

catch
    fprintf('encountered error, closing file\n');
    
end % try

fclose(fid);

end % ReadImageanalysisTxt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin)
% [Field, sParms] = z_GetDetectorViewer(zchan)
% [Field, sParms] = z_GetDetectorViewer(textfilename)
% [Field, sParms] = z_GetDetectorViewer(zchan, textfilename)
% [Field, sParms] = z_GetDetectorViewer(zchan, textfilename, settings.cfg)

% defaults:
zchan = [];
settingsfilename = [pwd '\ImaSettings.cfg'];
textfilename = [pwd '\imageanalysis.txt'];
settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename, see p. 634

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ltmp = SkipEmptyLines(fid)

    ltmp = [];
    while isempty(ltmp) && ~feof(fid),
        ltmp = fgetl(fid);
    end

end % SkipEmptyLines
