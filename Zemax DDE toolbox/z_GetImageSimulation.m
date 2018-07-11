function Sout = z_GetImageSimulation(varargin)
% Sout = z_GetImageSimulation(zchan)
% Sout = z_GetImageSimulation(textfilename)
% Sout = z_GetImageSimulation(zchan, textfilename)
% Sout = z_GetImageSimulation(zchan, textfilename, settings.cfg)
% 
% calls ZEMAX Image Simulation analysis
% the settings are the default settings for the lens
%
% Results are stored in ASCII form in textfilename, the
% default is [pwd '\imagesimu.txt']
% 
% efficiency = fraction of rays to reach the image plane
% img = [x,y] location of rays in the image plane
% field = [x,y] location of rays at source
% pupil = [x,y] location of rays at the pupil
% wave = wavelength# for each ray


[zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin{:});

if ~isempty(zchan),
    % Geometric Image Analysis
    typecode = 'Sim';
    
    % make the call
    cmdstr =...
        ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', ' settingsflag];
    retval = ddereq(zchan,cmdstr,[1 1]);
    
    % wait until zemax is ready
    while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
end

% get the image data from the textfile
Sout = ReadImageSimuTxt(textfilename);

% add x,y grid vectors for convenience
if isempty(Sout.ImageSize)
    % might be "Source Image"
    Sout.ImageSize = Sout.BitmapSize;
end
pixelsize = Sout.ImageSize ./ Sout.BitmapSize;
firstpixel = -Sout.ImageSize/2 + pixelsize/2;
lastpixel = -firstpixel;
Sout.x = linspace(firstpixel(1),lastpixel(1),Sout.BitmapSize(1))';
Sout.y = linspace(firstpixel(2),lastpixel(2),Sout.BitmapSize(2))';

end % main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sParms = ReadImageSimuTxt(textfile)

U = CConstants;

sParms = struct(...
    'ImageRGB',[]...
    ,'ImageIntensity',[]...
    ,'File', []...
    ,'Title', []...
    ,'Date', []...
    ,'Configuration', []...
    ,'DataType',[]...
    ,'AberrationsType',[]...
    ,'InputFile',[]...
    ,'BitmapSize',[]... % [height width] in pixels
    ,'ObjectHeight',[]...
    ,'FieldPosition',[]...
    ,'ImageSize',[]... % [height width] in (m)
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
                sParms.File = atmp{1}{1};
            case 'Title:',
                atmp = textscan(ltmp,'Title : %s','whitespace','\t');
                if ~isempty(atmp{1}), sParms.Title = atmp{1}{1}; end
            case 'Date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                sParms.Date = atmp{1}{1};
            case 'Configuration',
                atmp = textscan(ltmp,'Configuration %d');
                sParms.Configuration = atmp{1};
            case 'Data',
                sParms.DataType = ltmp(regexp(ltmp,': ')+1:end);
            case 'Aberrations',
                sParms.AberrationsType = ltmp(regexp(ltmp,': ')+1:end);
            case 'Input',
                sParms.InputFile = ltmp(regexp(ltmp,': ')+1:end);
            case 'Bitmap',
                % two lines, height and width
                htmp = sscanf(ltmp,'Bitmap Height : %d');
                ltmp = fgetl(fid);
                wtmp = sscanf(ltmp,'Bitmap Width  : %d');
                sParms.BitmapSize = [htmp wtmp];
            case 'Object',
                sParms.ObjectHeight = sscanf(ltmp,'Object Height : %f');
            case 'Field',
                sParms.FieldPosition = sscanf(ltmp,'Field position: %f');
            case 'Image',
                sParms.ImageSize = sscanf(ltmp(regexp(ltmp,':')+1:end),'%f W x %f H')*U.MM;
                sParms.ImageSize = sParms.ImageSize(:)'; % row vector
                
            case 'xpix',
                % data follows
                % format is:
                % xpix  ypix   R   G   B
                A = textscan(fid,'%d %d %d %d %d');
                [xpix, ypix, rr, gg, bb] = deal(A{:});
                sParms.ImageRGB = uint32(zeros(horzcat(sParms.BitmapSize, 3)));
                sParms.ImageRGB(:,:,1) = reshape(rr, sParms.BitmapSize);
                sParms.ImageRGB(:,:,2) = reshape(gg, sParms.BitmapSize);
                sParms.ImageRGB(:,:,3) = reshape(bb, sParms.BitmapSize);
                
                sParms.ImageIntensity = sum(sParms.ImageRGB,3);
            otherwise,

        end % switch

    end % while ~isempty
    
    % get next line
    ltmp = fgetl(fid);

end % while ~feof


catch
    fprintf('encountered error, closing file\n');
    
end % try

fclose(fid);

end % ReadImageSimuTxt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin)
% [Field, sParms] = z_GetDetectorViewer(zchan)
% [Field, sParms] = z_GetDetectorViewer(textfilename)
% [Field, sParms] = z_GetDetectorViewer(zchan, textfilename)
% [Field, sParms] = z_GetDetectorViewer(zchan, textfilename, settings.cfg)

% defaults:
zchan = [];
settingsfilename = [pwd '\ImaSettings.cfg'];
textfilename = [pwd '\imagesimu.txt'];
settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename, see p. 634

switch(nargin),
    case 0,
        [filename, pathname, fi] = uigetfile({'*.txt'},'Select Image Simulation Analysis');
        if ~isequal(filename,0),
            textfilename = [pathname filename];
        end

        % use defaults for everything else

    case 1,
        if ischar(varargin{1}),
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

