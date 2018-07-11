function Results = z_GetFiberCoupling(varargin)
% Results = z_GetFiberCoupling(zchan)
% Results = z_GetFiberCoupling(textfilename)
% Results = z_GetFiberCoupling(zchan, textfilename)
% Results = z_GetFiberCoupling
% 
% calls ZEMAX Fiber coupling analysis
% the settings are the default settings for the lens
%
% [...] = z_GetFiberCoupling(zchan)
%    ZEMAX calculates Fiber coupling on current lens, results are stored in
%    [pwd '\fclanalysis.txt']
%
% [...] = z_GetFiberCoupling(textfilename)
%    read results from textfilename
%
% [...] = z_GetFiberCoupling(zchan, textfilename)
%    ZEMAX calculates Fiber coupling on current lens and stores results in
%    textfilename
% 
% [...] = z_GetFiberCoupling
%    interactively select a text file with Fiber Coupling results to read
%
% Results = struct:
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
[zchan, textfilename] = ValidateInputs(varargin{:});

if ~isempty(zchan),
    % call ZEMAX to do the POP analysis
    
    % file for zemax settings
    %settingsfilename = [pwd '\FclSettings.cfg'];
    settingsfilename = ['c:' '\FclSettings.cfg'];
    % flag determines use/save settings. see discussion in "GetTextFile" in Zemax manual
    settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename

    % Analysis code
    typecode = 'Fcl';

    % make the call
    cmdstr =...
        ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', ' settingsflag];
    retval = ddereq(zchan,cmdstr,[1 1]);

    % wait until zemax is ready
    while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
end

% get the image data from the textfile
Results = ReadAnalysisTxt(textfilename);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function headerinfo = ReadAnalysisTxt(textfilename)

global NM UM MM;

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

% read each line
% not all header lines are always included, depending on settings
ltmp = fgetl(fid); % 
ltmp = fgetl(fid);

while ~feof(fid),
    ltmp = fgetl(fid);
    if ~isempty(ltmp),
        wordparse = textscan(ltmp,'%s');
    
        switch wordparse{1}{1}, % first word
            case 'File',
                atmp = textscan(ltmp,'File : %s','whitespace','\t');
                headerinfo.Lensfile = atmp{1}{1};
            case 'Title:',

            case 'Date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                headerinfo.Date = atmp{1}{1};
            case 'Wavelength',
                %headerinfo.Wavelength = sscanf(ltmp,'Wavelength : %f')*UM;
                headerinfo.Wavelength = str2num(wordparse{1}{3});
            case 'X',
                headerinfo.FieldX = str2num(wordparse{1}{5});
            case 'Y',
                headerinfo.FieldY = str2num(wordparse{1}{5});
            case 'Receiver',
                switch wordparse{1}{2}, % second word
                    case 'Efficiency',
                        headerinfo.ReceiverEff = str2num(wordparse{1}{4});
                    otherwise
                end
            case 'System',
                headerinfo.SystemEff = str2num(wordparse{1}{4});
            case 'Coupling',
                headerinfo.CouplingEff = str2num(wordparse{1}{4});
            case 'Maximum',
                headerinfo.MaximumEff = str2num(wordparse{1}{4});
            
            otherwise,

        end % switch
    end % if ~isempty line
end % while ~feof

fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zchan, textfilename] = ValidateInputs(varargin)
% Results = z_GetFiberCoupling(zchan)
% Results = z_GetFiberCoupling(textfilename)
% Results = z_GetFiberCoupling(zchan, textfilename)
% Results = z_GetFiberCoupling

switch(nargin),
    case 0,
        [filename, pathname, fi] = uigetfile({'*.txt'},'Select Fiber Coupling Analysis');
        if isequal(filename,0),
            textfilename = [];
        else,
            textfilename = [pathname filename];
        end
        zchan = [];
    case 1,
        if isstr(varargin{1}),
            textfilename = varargin{1};
            zchan = [];
        else
            zchan = varargin{1};
            textfilename = [pwd '\fclanalysis.txt'];
        end
    case 2,
        zchan = varargin{1};
        textfilename = varargin{2};
    otherwise,
        error('usage: [Field, parmsstruct] = z_GetFiberCoupling(zchan, textfilename)');
end

end % ValidateInputs
