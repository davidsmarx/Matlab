function [r, psf, headerinfo] = z_GetPSFFFTcrosssection(varargin)
% [r, PSF, parmsstruct] = z_GetPSFFFTcrosssection(zchan)
% [r, PSF, parmsstruct] = z_GetPSFFFTcrosssection(textfilename)
% [r, PSF, parmsstruct] = z_GetPSFFFTcrosssection(zchan, textfilename)
% [r, PSF, parmsstruct] = z_GetPSFFFTcrosssection
% 
% calls ZEMAX FFT PSF analysis
% the settings are the default settings for the lens
%
% [...] = z_GetPSFFFTcrosssection(zchan)
%    ZEMAX calculates PSF on current lens, results are stored in
%    [pwd '\PSFFFTanalysis.txt']
%
% [...] = z_GetPSFFFTcrosssection(textfilename)
%    read results from textfilename
%
% [...] = z_GetPSFFFTcrosssection(zchan, textfilename)
%    ZEMAX calculates POP on current lens and stores results in
%    textfilename
% 
% [...] = z_GetPSFFFTcrosssection
%    interactively select a text file with POP results to read
%
% psf = vector
% r   = radius (x)
% headerinfo:
%    Lensfile
%    Date
%    Datatype (huygens psf, etc.)
%    Field (position mm)
%    Wavelength = [start_wave end_wave]
%    NxNyPupil = [Nx, Ny] grid size at pupil
%    NxNyImage = [Nx, Ny] grid size at image

% validate inputs
[zchan, textfilename] = ValidateInputs(varargin{:});
if isempty(textfilename), [r, psf, headerinfo] = deal([]); return, end

if ~isempty(zchan),
    % call ZEMAX to do the analysis
    
    % zemax three-letter code for this analysis
    typecode = 'Pcs';
    
    % file for zemax settings
    settingsfilename = [pwd '\' typecode 'Settings.cfg'];
    % flag determines use/save settings. see discussion in "GetTextFile" in Zemax manual
    settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename

    % make the call
    cmdstr =...
        ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', ' settingsflag];
    retval = ddereq(zchan,cmdstr,[1 1]);

    % wait until zemax is ready
    while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
end

% get the image data from the textfile
[r, psf, headerinfo] = ReadAnalysisTxt(textfilename);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, psf, headerinfo] = ReadAnalysisTxt(textfilename)

global NM UM MM;

fid = fopen(textfilename,'rt');
if fid == -1, error(['error opening text file ' textfilename]); end

% psf = vector
% r   = radius (x)
% headerinfo:
%    Lensfile
%    Date
%    Datatype (huygens psf, etc.)
%    Field (position)
%    Wavelength = [start_wave end_wave]
%    NxNyPupil = [Nx, Ny] grid size at pupil
%    NxNyImage = [Nx, Ny] grid size at image

% read each line, assumes 'Values...' is the last header line before data
% not all header lines are always included, depending on settings
ltmp = fgetl(fid); % 
ltmp = fgetl(fid);

while ~feof(fid),
    ltmp = fgetl(fid);
    if ~isempty(ltmp) & (ltmp ~= -1),
        wordparse = textscan(ltmp,'%s');
    
        switch lower(wordparse{1}{1}), % first word
            case 'file',
                atmp = textscan(ltmp,'File : %s','whitespace','\t');
                headerinfo.Lensfile = atmp{1}{1};
            case 'title:',

            case 'date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                headerinfo.Date = atmp{1}{1};
                
            case 'data',
                switch lower(wordparse{1}{2}),
                    case 'type',
                        headerinfo.Datatype = sscanf(ltmp,'Data Type : %s');
                    case 'position',
                        % read the data
                        ltmp = fgetl(fid);
                        atmp = fscanf(fid,'%f %f %f',[3 inf])';

                        [r, psf] = deal(atmp(:,2)*UM,atmp(:,3));
                end
                
            case 'field',
                headerinfo.Field = sscanf(ltmp,'Field : %f')*MM;
                
            case 'wavelength',
                headerinfo.Wavelength = sscanf(ltmp,'Wavelength : %f to %f')*UM;
                
            case 'pupil',
                headerinfo.NxNyPupil = sscanf(ltmp,'Pupil grid size: %d by %d');

            case 'image',
                headerinfo.NxNyImage = sscanf(ltmp,'Image grid size: %d by %d');

            case 'values',
                % now read data
                
        end % switch
    end % if ~isempty line
end % while ~feof

fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zchan, textfilename] = ValidateInputs(varargin)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan)
% [Field, parmsstruct] = z_GetPOPAnalysis(textfilename)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan, textfilename)

switch(nargin),
    case 0,
        [filename, pathname, fi] = uigetfile({'*.txt'},'Select POP Analysis');
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
            textfilename = [pwd '\popanalysis.txt'];
        end
    case 2,
        zchan = varargin{1};
        textfilename = varargin{2};
    otherwise,
        error('usage: [Field, parmsstruct] = z_GetPOPAnalysis(zchan, textfilename)');
end

end % ValidateInputs
