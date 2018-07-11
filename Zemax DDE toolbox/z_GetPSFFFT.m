function [Field, headerinfo] = z_GetPSFFFT(varargin)
% [PSF, parmsstruct] = z_GetPSFFFT(zchan)
% [PSF, parmsstruct] = z_GetPSFFFT(textfilename)
% [PSF, parmsstruct] = z_GetPSFFFT(zchan, textfilename)
% [PSF, parmsstruct] = z_GetPSFFFT
% 
% calls ZEMAX FFT PSF analysis
% the settings are the default settings for the lens
%
% [...] = z_GetPSFFFT(zchan)
%    ZEMAX calculates PSF on current lens, results are stored in
%    [pwd '\PSFFFTanalysis.txt']
%
% [...] = z_GetPSFFFT(textfilename)
%    read results from textfilename
%
% [...] = z_GetPSFFFT(zchan, textfilename)
%    ZEMAX calculates POP on current lens and stores results in
%    textfilename
% 
% [...] = z_GetPSFFFT
%    interactively select a text file with POP results to read
%
% PSF = matrix of PSF values
% headerinfo:
%    Lensfile
%    Date
%    Datatype (polychromatic fft psf, etc.)
%    Wavelength = [start_wave end_wave]
%    Pointspacing = [dx dy]
%    Dataarea
%    Surface
%    ReferenceCoord = [x0, y0]
%    NxNyPupil = [Nx, Ny] grid size at pupil
%    NxNyImage = [Nx, Ny] grid size at image
%    Center.ir = row index of center point
%    Center.ic = column index of center point

% validate inputs
[zchan, textfilename] = ValidateInputs(varargin{:});

if ~isempty(zchan),
    % call ZEMAX to do the analysis
    
    % file for zemax settings
    settingsfilename = [pwd '\FpsSettings.cfg'];
    % flag determines use/save settings. see discussion in "GetTextFile" in Zemax manual
    settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename

    % Analysis code
    typecode = 'Fps';

    % make the call
    cmdstr =...
        ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', ' settingsflag];
    retval = ddereq(zchan,cmdstr,[1 1]);

    % wait until zemax is ready
    while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
end

% get the image data from the textfile
[Field, headerinfo] = ReadAnalysisTxt(textfilename);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Field, headerinfo] = ReadAnalysisTxt(textfilename)

global NM UM MM;

fid = fopen(textfilename,'rt');
if fid == -1, error(['error opening text file ' textfilename]); end

% read header info
%    Lensfile
%    Date
%    Datatype (polychromatic fft psf, etc.)
%    Wavelength = [start_wave end_wave]
%    Pointspacing = [dx dy]
%    Dataarea
%    Surface
%    ReferenceCoord = [x0, y0]
%    NxNyPupil = [Nx, Ny] grid size at pupil
%    NxNyImage = [Nx, Ny] grid size at image
%    Center.ir = row index of center point
%    Center.ic = column index of center point

% read each line, assumes 'Values...' is the last header line before data
% not all header lines are always included, depending on settings
ltmp = fgetl(fid); % 
ltmp = fgetl(fid);

while ~feof(fid),
    if ~isempty(ltmp)
        wordparse = textscan(ltmp,'%s');
    
        switch lower(wordparse{1}{1}), % first word
            case 'file',
                atmp = textscan(ltmp,'File : %s','whitespace','\t');
                headerinfo.Lensfile = atmp{1}{1};
            case 'title:',

            case 'date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                headerinfo.Date = atmp{1}{1};
            case 'polychromatic',
                headerinfo.Datatype = ltmp;
            case 'data',
                switch lower(wordparse{1}{2})
                    case 'spacing',
                        atmp = sscanf(ltmp,'Data spacing is %f');
                        headerinfo.Pointspacing = [atmp atmp]*UM;
                    case 'area',
                        atmp = sscanf(ltmp,'Data area is %f');
                        headerinfo.Dataarea = atmp*UM;
                end
            case 'surface',
                atmp = sscanf(ltmp,'Surface: %d');
                if isempty(atmp),
                    % must be Surface: Image
                    headerinfo.Surface = inf;
                else
                    headerinfo.Surface = atmp;
                end
            case 'reference',
                headerinfo.ReferenceCoord = sscanf(ltmp,'Reference Coordinates: %f, %f');
                
            case 'pupil',
                headerinfo.NxNyPupil = sscanf(ltmp,'Pupil grid size: %d by %d');
            case 'image',
                headerinfo.NxNyImage = sscanf(ltmp,'Image grid size: %d by %d');
                
            case 'center',
                atmp = sscanf(ltmp,'Center point is: row %d, column %d');
                headerinfo.Center.ir = atmp(1);
                headerinfo.Center.ic = atmp(2);

            case 'values',
                % now read data
                ltmp = fgetl(fid); % blank line before data
                % read the data
                Field = fscanf(fid,'%e',headerinfo.NxNyImage')';
                
            otherwise,
                % header lines that begin with the number value must be
                % deciphered by the word at the end of the line
                switch lower(wordparse{1}{end}),
                    case 'deg.',
                        atmp = sscanf(ltmp,'%f to %f');
                        headerinfo.Wavelength = atmp;
                end

        end % switch
    end % if ~isempty line

    % read next line
    ltmp = fgetl(fid);

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
