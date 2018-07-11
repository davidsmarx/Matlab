function [Field, headerinfo] = z_GetPOPAnalysis(varargin)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan)
% [Field, parmsstruct] = z_GetPOPAnalysis(textfilename)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan, textfilename)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan, textfilename, settings.cfg)
% [Field, parmsstruct] = z_GetPOPAnalysis
% 
% calls ZEMAX physical optics propagation
% the settings are the default settings for the lens
%
% [...] = z_GetPOPAnalysis(zchan)
%    ZEMAX calculates POP on current lens, results are stored in
%    [pwd '\popnalysis.txt']
%
% [...] = z_GetPOPAnalysis(textfilename)
%    read results from textfilename
%
% [...] = z_GetPOPAnalysis(zchan, textfilename)
%    ZEMAX calculates POP on current lens and stores results in
%    textfilename
% 
% [...] = z_GetPOPAnalysis
%    interactively select a text file with POP results to read
%
% [...] = z_GetPOPAnalysis(zchan, textfilename, settings.cfg)
%    use the configuration in settings.cfg to make the calculations. If
%    settings.cfg is not specified, default settings are used, and then the
%    default settings are written to 'popsettings.cfg'
%
% Field = matrix of field values calculated from the POP
% headerinfo:
%    Lensfile
%    Date
%    Datatype (Irradience, etc.)
%    NxNy = [Nx, Ny] array size
%    Pointspacing = [dx dy]
%    Wavelength
%    PilotSize
%    PilotWaist
%    PilotPos
%    PilotRayleigh

% validate inputs
[zchan, textfilename, settingsfilename, settingsflag] =...
    ValidateInputs(varargin{:});

if ~isempty(zchan),
    % call ZEMAX to do the POP analysis
    
    % Analysis code
    typecode = 'Pop';
    
    % make the call
    cmdstr =...
        ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', ' settingsflag];
    retval = ddereq(zchan,cmdstr,[1 1]);

    % wait until zemax is ready
    while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
end

if ~isempty(textfilename),
    % get the image data from the textfile
    [Field, headerinfo] = ReadAnalysisTxt(textfilename);
    
    % convert Field units to MKS
    Field = CheckFieldUnits(Field, headerinfo);
    
    % add grid arrays
    headerinfo = AppendPixelGrid(headerinfo);
    
else,
    % maybe zbf file is used
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Field = CheckFieldUnits(Fieldtmp, headerinfo)

U = CConstants;

% convert units of Field to MKS
if ~isempty(strfind(headerinfo.Datatype, 'Total Irradiance')),
    Field = Fieldtmp/U.MM/U.MM; % convert W/mm^2 to W/m^2
else
    Field = Fieldtmp;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function headerinfo = AppendPixelGrid(htmp)

headerinfo = htmp;
ix = -htmp.NxNy(1)/2:(htmp.NxNy(1)/2-1);
iy = -htmp.NxNy(2)/2:(htmp.NxNy(2)/2-1);

headerinfo.xg = headerinfo.Pointspacing(1) * ix(:);
headerinfo.yg = headerinfo.Pointspacing(2) * iy(:);

end % AppendPixelGrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Field, headerinfo] = ReadAnalysisTxt(textfilename)

U = CConstants;

fid = fopen(textfilename,'rt');

if fid == -1, error(['error opening text file ' textfilename]); end

% read header info
% headerinfo:
%    Lensfile
%    Date
%    Datatype (Irradience, etc.)
%    Pointspacing = [dx dy]
%    Wavelength
%    FiberCoupling: [system, receiver, coupling]
%    Pilot: [Size PilotWaist PilotPos PilotRayleigh]

% read each line, assumes 'Pilot' is the last header line before data
% not all header lines are always included, depending on settings
ltmp = fgetl(fid);

while ~feof(fid)
    if ~isempty(ltmp)
        wordparse = textscan(ltmp,'%s');

        switch wordparse{1}{1}, % first word
            case 'File',
                atmp = textscan(ltmp,'File : %s','whitespace','\t');
                headerinfo.Lensfile = atmp{1}{1};
            case 'Title:',

            case 'Date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                headerinfo.Date = atmp{1}{1};
            case {'Total', 'Phase'}
                atmp = textscan(ltmp,'%s[^\n]','whitespace','\t');
                headerinfo.Datatype = atmp{1}{1};
            case 'Grid',
                headerinfo.NxNy = sscanf(ltmp,'Grid size (X by Y): %d by %d');
            case 'Point',
                headerinfo.Pointspacing = sscanf(ltmp,'Point spacing (X by Y): %f by %f')'*U.MM;
            case 'Wavelength',
                % Wavelength 0.66466 µm in index 1.00000 at 0.0000, -0.0250 mm
                % apply units
                A = textscan(ltmp,'Wavelength %f %s in index %f at %f, %f %s');
                headerinfo.Wavelength = A{1}*U.UM;
                headerinfo.index = A{3};
                headerinfo.Field = [A{4} A{5}]*U.MM;
            case 'Display',
                % total grid size = nx * pointspacing
                atmp = textscan(ltmp,'Display X Width = %f, Y Height = %f'); % millimeters
                headerinfo.DisplaySize = [
                    atmp{1}*U.MM atmp{2}*U.MM
                    ]; % x,y

            case 'Peak',
                atmp = sscanf(ltmp,'Peak Irradiance = %f Watts/Millimeters^2, Total Power = %f Watts'); % peak irradiance
                headerinfo.PeakIrradiance = atmp(1)/U.MM/U.MM;
                headerinfo.TotalPower = atmp(2); % Watts
            case 'Fiber',
                atmp = sscanf(ltmp,'Fiber Efficiency: System %f, Receiver %f, Coupling %f');
                headerinfo.FiberCoupling.system = atmp(1);
                headerinfo.FiberCoupling.receiver = atmp(2);
                headerinfo.FiberCoupling.total = atmp(3);
                
            case 'Pilot:',
                headerinfo.Pilot = sscanf(ltmp,'Pilot: Size= %f, Waist= %f, Pos= %f, Rayleigh= %f')';

                % now read data (assumes Pilot is the last header line before
                % data)
                ltmp = fgetl(fid); % blank line before data
                % read the data
                Field = fscanf(fid,'%e',headerinfo.NxNy')';
            otherwise,

        end % switch

    end % while ~isempty
    
    % get next line
    ltmp = fgetl(fid);

end % while ~feof

fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan)
% [Field, parmsstruct] = z_GetPOPAnalysis(textfilename)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan, textfilename)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan, textfilename, settings.cfg)

switch(nargin),
    case 0,
        [filename, pathname, fi] = uigetfile({'*.txt'},'Select POP Analysis');
        if isequal(filename,0),
            textfilename = [];
        else,
            textfilename = [pathname filename];
        end
        zchan = [];
        settingsfilename = [pwd '\PopSettings.cfg'];

        % flag determines use/save settings. see discussion in "GetTextFile" in Zemax manual
        settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename

    case 1,
        if isstr(varargin{1}),
            textfilename = varargin{1};
            zchan = [];
        else
            zchan = varargin{1};
            textfilename = [pwd '\popanalysis.txt'];
        end
        settingsfilename = [pwd '\PopSettings.cfg'];
        % flag determines use/save settings. see discussion in "GetTextFile" in Zemax manual
        settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename
    case 2,
        zchan = varargin{1};
        textfilename = varargin{2};
        settingsfilename = [pwd '\PopSettings.cfg'];
        % flag determines use/save settings. see discussion in "GetTextFile" in Zemax manual
        settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename
    case 3,
        [zchan, textfilename, settingsfilename] = deal(varargin{:});
        settingsflag = '1'; % 1 => use settings in settingsfilename
    otherwise,
        error('usage: [Field, parmsstruct] = z_GetPOPAnalysis(zchan, textfilename)');
end

end % ValidateInputs
