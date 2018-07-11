function [Wavefront, sParms] = z_GetWavefrontMap(varargin)
% [Wavefront, sParms] = z_GetWavefrontMap(zchan, textfilename)
% 
% [Wavefront, sParms] = z_GetWavefrontMap(zchan)
% [Wavefront, sParms] = z_GetWavefrontMap(textfilename)
% [Wavefront, sParms] = z_GetWavefrontMap(zchan, textfilename)
% [Wavefront, sParms] = z_GetWavefrontMap(zchan, textfilename, settings.cfg)
% [Wavefront, sParms] = z_GetWavefrontMap
%
% calls ZEMAX Wavefront -> Wavefrontmap analysis
% the settings are the default settings for the lens
%
% [...] = z_GetWavefrontMap(zchan)
%    ZEMAX calculates the Wavefrontmap on current lens, results are stored in
%    [pwd '\wavefrontmap.txt']
%
% [...] = z_GetWavefrontMap(textfilename)
%    read results from textfilename
%
% [...] = z_GetWavefrontMap(zchan, textfilename)
%    ZEMAX calculates the wavefrontmap on current lens and stores results in
%    textfilename
% 
% [...] = z_GetWavefrontMap
%    interactively select a text file with wavefrontmap results to read
%
% [...] = z_GetWavefrontMap(zchan, textfilename, settings.cfg)
%    use the configuration in settings.cfg to make the calculations. If
%    settings.cfg is not specified, default settings are used, and then the
%    default settings are written to 'popsettings.cfg'
%
% Wavefront = Wavefront
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

if ~isempty(zchan),
    % Geometric Image Analysis
    typecode = 'Wfm';
    
    % make the call
    cmdstr =...
        ['GetTextFile, "' textfilename '", ' typecode ', "' settingsfilename '", ' settingsflag];
    retval = ddereq(zchan,cmdstr,[1 1]);
    
    % wait until zemax is ready
    while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
end

% get the image data from the textfile
[Wavefront, sParms] = ReadWavefrontmap(textfilename);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Wavefront, sParms] = ReadWavefrontmap(textfile)

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
Wavefront = [];
sParms = struct;

fid = fopen(textfile,'rt','native','csUnicode');
if fid == -1, error(['error opening text file ' textfile]); end

% try,
    
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
                if ~isempty(atmp{1}),
                    sParms.Title = atmp{1}{1};
                end
            case 'Date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                sParms.Date = atmp{1}{1};
            case 'Wavefront',
                %Wavefront Function
                %0.5448 µm at 0.0000, 0.0000 mm.

                % next line is wavelength and field point
                ltmp = fgetl(fid);
                while isempty(ltmp), ltmp = fgetl(fid); end
                % apply units
                A = textscan(ltmp,'%f %s at %f, %f %s');
                sParms.Wavelength = A{1}*U.UM;
                sParms.Field = [A{3} A{4}]*U.MM;

            case 'Peak',
                % Peak to Valley = 1.1354 waves, Fringes/Wave = 1.0000.
                A = textscan(ltmp,'Peak to valley = %f waves, RMS = %f');
                sParms.PeakToValley = A{1};
                sParms.RMS = A{2};
                                
            case 'Surface:',
                % Surface: Image
                A = textscan(ltmp,'Surface: %s');
                sParms.Surface = A{1}{1};
                
            case 'Exit'
                % Exit Pupil Diameter: 7.6531E+000 Millimeters
                A = textscan(ltmp, 'Exit Pupil Diameter: %f');
                sParms.ExitPupilDiameter = A{1}*U.MM;
                
            case 'Pupil'
                % Pupil grid size:
                A = textscan(ltmp, 'Pupil grid size: %f by %f');
                sParms.GridSize = [A{1} A{2}]; % [x, y] pixels
                
            case 'Center'
                % Center point is:
                A = textscan(ltmp, 'Center point is: Col %f, Row %f');
                sParms.GridCenter = [A{1} A{2}]; % [x0, y0] pixels
                
                % data follows
                ltmp = fgetl(fid); % empty line
                Data = fscanf(fid,'%f');
                Wavefront = reshape(Data,sParms.GridSize(2),sParms.GridSize(1));
                
            otherwise,

        end % switch

    end % while ~isempty
    
    % get next line
    ltmp = fgetl(fid);

end % while ~feof

% catch,
%     fprintf('encountered error, closing file\n');
%     
% end

fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin)
% [Field, sParms] = z_GetWavefrontMap(zchan)
% [Field, sParms] = z_GetWavefrontMap(textfilename)
% [Field, sParms] = z_GetWavefrontMap(zchan, textfilename)
% [Field, sParms] = z_GetWavefrontMap(zchan, textfilename, settings.cfg)

switch(nargin),
    case 0,
        [filename, pathname, fi] = uigetfile({'*.txt'},'Select wavefrontmap Analysis');
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
        error('usage: [Field, sParms] = z_GetWavefrontMap(zchan, textfilename)');
end

end % ValidateInputs
