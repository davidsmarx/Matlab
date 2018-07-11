function [fringeamp, sParms] = z_GetInterferogram(varargin)
% [fringeamp, sParms] = z_GetInterferogram(zchan, textfilename)
% 
% [Field, sParms] = z_GetInterferogram(zchan)
% [Field, sParms] = z_GetInterferogram(textfilename)
% [Field, sParms] = z_GetInterferogram(zchan, textfilename)
% [Field, sParms] = z_GetInterferogram(zchan, textfilename, settings.cfg)
% [Field, sParms] = z_GetInterferogram
%
% calls ZEMAX Wavefront -> Interferogram analysis
% the settings are the default settings for the lens
%
% [...] = z_GetInterferogram(zchan)
%    ZEMAX calculates the Interferogram on current lens, results are stored in
%    [pwd '\interferogram.txt']
%
% [...] = z_GetInterferogram(textfilename)
%    read results from textfilename
%
% [...] = z_GetInterferogram(zchan, textfilename)
%    ZEMAX calculates the Interferogram on current lens and stores results in
%    textfilename
% 
% [...] = z_GetInterferogram
%    interactively select a text file with Interferogram results to read
%
% [...] = z_GetInterferogram(zchan, textfilename, settings.cfg)
%    use the configuration in settings.cfg to make the calculations. If
%    settings.cfg is not specified, default settings are used, and then the
%    default settings are written to 'popsettings.cfg'
%
% fringeamp = fringe amplitude 0 to 1.0
% sParms fields:
%    file
%    date
%    configuration
%    description
%    wavelength
%    field
%    surface
%    exitpupildiam
%    xtilt
%    ytilt

% validate inputs
[zchan, textfilename, settingsfilename, settingsflag] =...
    ValidateInputs(varargin{:});

% Geometric Image Analysis
typecode = 'Int';

% make the call
cmdstr =...
    ['GetTextFile, "' textfilename '", ' typecode ', "' settingsfilename '", ' settingsflag];
retval = ddereq(zchan,cmdstr,[1 1]);

% wait until zemax is ready
while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end

% get the image data from the textfile
[fringeamp, sParms] = ReadInterferogram(textfilename);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fringeamp, sParms] = ReadInterferogram(textfile)

U = CConstants;

% sParms fields:
%    file
%    date
%    configuration
%    description
%    wavelength
%    field
%    surface
%    exitpupildiam
%    xtilt
%    ytilt

% initialize outputs
fringeamp = [];
sParms = struct;

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

            case 'Date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                sParms.Date = atmp{1}{1};
            case 'Interferogram',
                %Interferogram between Configurations 1 and 2
                %0.5448 µm at 0.0000, 0.0000 mm.

                % next line is wavelength and field point
                ltmp = fgetl(fid);
                % apply units
                A = textscan(ltmp,'%f %s at %f, %f %s');
                sParms.Wavelength = A{1}*U.UM;
                sParms.Field = [A{3} A{4}]*U.MM;

            case 'Peak',
                % Peak to Valley = 1.1354 waves, Fringes/Wave = 1.0000.
                A = textscan(ltmp,'Peak to Valley = %f waves, Fringes/Wave = %f');
                sParms.PeakToValley = A{1};
                sParms.FringesPerWave = A{2};
                                
            case 'Surface:',
                % Surface: Image
                A = textscan(ltmp,'Surface: %s');
                sParms.Surface = A{1}{1};
                
            case 'Exit'
                % Exit Pupil Diameter: 7.6531E+000 Millimeters
                A = textscan(ltmp, 'Exit Pupil Diameter: %f');
                sParms.ExitPupilDiameter = A{1}*U.MM;
                
            case 'Xtilt'
                % Xtilt = 0.00, Ytilt = 0.00.
                A = textscan(ltmp, 'Xtilt = %f, Ytilt = %f');
                sParms.XYtilt = [A{1} A{2}];
                
                % data follows
                ltmp = fgetl(fid); % empty line
                Data = fscanf(fid,'%f');
                M = length(Data);
                N = sqrt(M);
                fringeamp = reshape(Data,N,N);
                
            otherwise,

        end % switch

    end % while ~isempty
    
    % get next line
    ltmp = fgetl(fid);

end % while ~feof

catch,
    fprintf('encountered error, closing file\n');
    
end

fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin)
% [Field, sParms] = z_GetInterferogram(zchan)
% [Field, sParms] = z_GetInterferogram(textfilename)
% [Field, sParms] = z_GetInterferogram(zchan, textfilename)
% [Field, sParms] = z_GetInterferogram(zchan, textfilename, settings.cfg)

switch(nargin),
    case 0,
        [filename, pathname, fi] = uigetfile({'*.txt'},'Select Interferogram Analysis');
        if isequal(filename,0),
            textfilename = [];
        else,
            textfilename = [pathname filename];
        end
        zchan = [];
        settingsfilename = [pwd '\IntSettings.cfg'];

        % flag determines use/save settings. see discussion in "GetTextFile" in Zemax manual
        settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename

    case 1,
        if isstr(varargin{1}),
            textfilename = varargin{1};
            zchan = [];
        else
            zchan = varargin{1};
            textfilename = [pwd '\intanalysis.txt'];
        end
        settingsfilename = [pwd '\IntSettings.cfg'];
        % flag determines use/save settings. see discussion in "GetTextFile" in Zemax manual
        settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename
    case 2,
        zchan = varargin{1};
        textfilename = varargin{2};
        settingsfilename = [pwd '\IntSettings.cfg'];
        % flag determines use/save settings. see discussion in "GetTextFile" in Zemax manual
        settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename
    case 3,
        [zchan, textfilename, settingsfilename] = deal(varargin{:});
        settingsflag = '1'; % 1 => use settings in settingsfilename
    otherwise,
        error('usage: [Field, sParms] = z_GetInterferogram(zchan, textfilename)');
end

end % ValidateInputs
