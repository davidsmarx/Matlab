function [Field, headerinfo] = z_GetPSFHuygens(varargin)
% [PSF, parmsstruct] = z_GetPSFHuygens(zchan, options)
% [PSF, parmsstruct] = z_GetPSFHuygens(options)
% 
% options:
%     'textfilename', textfilename
%     'settingsfilename', settingsfilename ( .cfg )
%
% calls ZEMAX FFT PSF analysis
% the settings are the default settings for the lens
%
% [...] = z_GetPSFHuygens(zchan)
%    ZEMAX calculates PSF on current lens, results are stored in
%    [pwd '\PSFFFTanalysis.txt']
%
% [...] = z_GetPSFHuygens(textfilename)
%    read results from textfilename
%
% [...] = z_GetPSFHuygens(zchan, textfilename)
%    ZEMAX calculates POP on current lens and stores results in
%    textfilename
% 
% [...] = z_GetPSFHuygens
%    interactively select a text file with POP results to read
%
% Field = matrix of field values calculated from the POP
% headerinfo:
%    Lensfile
%    Date
%    Datatype (huygens psf, etc.)
%    Wavelength = [start_wave end_wave]
%    Pointspacing = [dx dy]
%    Dataarea
%    StrehlRatio
%    NxNyPupil = [Nx, Ny] grid size at pupil
%    NxNyImage = [Nx, Ny] grid size at image
%    Center.ir = row index of center point
%    Center.ic = column index of center point
%    ReferenceCoord = [x0, y0]

% validate inputs
[zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin{:});

if ~isempty(zchan),
    % call ZEMAX to do the analysis
    
    % zemax three-letter code for this analysis
    typecode = 'Hps';
    
    % make the call
    cmdstr =...
        ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', ' settingsflag];
    retval = ddereq(zchan,cmdstr,[1 1]);

    % wait until zemax is ready
    while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end
end

% get the image data from the textfile
[Field, headerinfo] = ReadAnalysisTxt(textfilename);
headerinfo.textfilename = textfilename;
headerinfo.settingsfilename = settingsfilename;
headerinfo.settingsflag = settingsflag;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Field, headerinfo] = ReadAnalysisTxt(textfilename)

global NM UM MM;

fid = fopen(textfilename,'rt');
if fid == -1, error(['error opening text file ' textfilename]); end

% read header info
% headerinfo:
%    Lensfile
%    Date
%    Datatype (huygens psf, etc.)
%    Wavelength = [start_wave end_wave]
%    Pointspacing = [dx dy]
%    Dataarea
%    StrehlRatio
%    NxNyPupil = [Nx, Ny] grid size at pupil
%    NxNyImage = [Nx, Ny] grid size at image
%    Center.ir = row index of center point
%    Center.ic = column index of center point
%    ReferenceCoord = [x0, y0]

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

            case 'date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                headerinfo.Date = atmp{1}{1};

            case 'huygens',
                headerinfo.Datatype = ltmp;

            case 'configuration',
                headerinfo.Configuraton = str2double(wordparse{1}([2 4]));
                % a couple of blank lines, then wavelength and field point
                ltmp = fgetl(fid);
                ltmp = fgetl(fid);
                ltmp = fgetl(fid);
                ctmp = textscan(ltmp,'%s');
                headerinfo.wave = str2double(ctmp{1}{1})*UM;
                headerinfo.fieldpoint = str2double(ctmp{1}(4:5))*MM;
                
            case 'data',
                switch lower(wordparse{1}{2})
                    case 'spacing',
                        atmp = sscanf(ltmp,'Data spacing is %f');
                        headerinfo.Pointspacing = [atmp atmp]*UM;
                    case 'area',
                        atmp = sscanf(ltmp,'Data area is %f');
                        headerinfo.Dataarea = atmp*UM;
                end
                
                
            case 'strehl',
                headerinfo.StrehlRatio = sscanf(ltmp,'Strehl ratio: %f');
                
            case 'pupil',
                headerinfo.NxNyPupil = sscanf(ltmp,'Pupil grid size: %d by %d')';
            case 'image',
                headerinfo.NxNyImage = sscanf(ltmp,'Image grid size: %d by %d')';

            case 'center',
                switch lower(wordparse{1}{2}),
                    case 'point'
                        atmp = sscanf(ltmp,'Center point is: %d, %d');
                        headerinfo.Center.ir = atmp(1);
                        headerinfo.Center.ic = atmp(2);
                    case 'coordinates'
                        headerinfo.CenterCoords = str2double(wordparse{1}(4:5));

                end

            case 'centroid'
                switch lower(wordparse{1}{2})
                    case 'offset'
                        headerinfo.CentroidOffset = str2double(wordparse{1}(4:5));
                    case 'coordinates:'
                        headerinfo.CentroidCoords = str2double(wordparse{1}(3:4));
                end
                
            case 'values',
                % now read data
                headerinfo.Values = ltmp;
                ltmp = fgetl(fid); % blank line before data
                % read the data
                Field = fscanf(fid,'%e',headerinfo.NxNyImage)';
                
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
function [zchan, textfilename, settingsfilename, settingsflag] = ValidateInputs(varargin)
% [Field, parmsstruct] = z_GetPOPAnalysis(zchan, options)
% [Field, parmsstruct] = z_GetPOPAnalysis(options)

zchan = [];
if isa(varargin{1},'double'),
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
