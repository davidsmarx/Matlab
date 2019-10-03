function [map, sParms] = ReadCodeVINTdata(fn, varargin)
% [map, sParms] = ReadCodeVINTdata(fn, varargin)
%
% map is wavefront in [radians]
%
% options:
%    'plot' (false) or true

U = CConstants;

% options
bPlot = CheckOption('plot', false, varargin{:});

% open
fid = fopen(fn,'rt');
if isequal(fid,-1),
    error(['cannot open file, ' fn]);
end

sParms = sParmsBlank;

% first valid line is title/description
sParms.Title = GetNextLine(fid); % skips comments and blank lines

% second line is Parameters
sParms.strParms = GetNextLine(fid);

% parse parameters    
wtmp = textscan(sParms.strParms,'%s');
wordparse = wtmp{1};
while ~isempty(wordparse),
    switch lower(wordparse{1}),
        case 'grd'
            sParms.DataType = 'GridData';
            sParms.ncols = str2double(wordparse{2});
            sParms.nrows = str2double(wordparse{3});
            wordparse = wordparse(4:end);
            
        case 'zfr' % zernike fringe coefficients
            sParms.DataType = 'ZernikeFringe';
            sParms.nz = str2double(wordparse{2});
            sParms.nrows = sParms.nz;
            sParms.ncols = 1;
            wordparse = wordparse(3:end);
            
        case 'wfr'
            sParms.DataValue = 'wavefront';
            wordparse = wordparse(2:end);
            
        case 'wvl' % wavelength in microns
            sParms.wavelength = str2double(wordparse{2})*U.UM;
            wordparse = wordparse(3:end);
            
        case 'ssz' % data value equal to one wave
            sParms.datascale_wave = str2double(wordparse{2});
            wordparse = wordparse(3:end);
            
        case 'nda' % no data value
            sParms.nda = str2double(wordparse{2});
            wordparse = wordparse(3:end);
            
        otherwise
            wordparse = wordparse(2:end);
            
    end % switch wordparse
    
end % while ~isempty wordparse

% data follows
% from Code V Help "Using Interferometric Data":
%   Each value is entered as an integer and it is assumed to be read by rows
%   from -X to +X (left to right) starting at +Y (top). 
data = fscanf(fid,'%f',[sParms.ncols sParms.nrows]);
fclose(fid);

if strcmp(sParms.DataType, 'GridData'),
    data = flipud(data.'); %
    data(data==sParms.nda) = nan;
end
if strcmp(sParms.DataType, 'ZernikeFringe'),
    % normalize zernike coeff to rms
    %rmsnorm = [1 2 2 sqrt(3) sqrt(6)*[1 1] sqrt(8)*[1 1 1 1] sqrt(5) sqrt(10)*[1 1 1 1] sqrt(3)*ones(1,6) sqrt(7) sqrt(14)*ones(1,6) 4*ones(1,8)];
    rmsnorm  =1;
    data = data./rmsnorm;
end

% calibrate, return radians
%map = 2*pi*data/sParms.datascale_wave; % return map in radians
% return waves, Code V default
map = data/sParms.datascale_wave;

if bPlot,
    figure, imageschcit(map)
    colorbartitle('Wavefront (rad)')
    title([pwd2titlestr(sParms.Title) '; \lambda = ' num2str(sParms.wavelength/U.NM) 'nm'])
end % if plot

end % main

function tline = GetNextLine(fid)

tline = '';

while ~feof(fid),
    tline = fgetl(fid);
    if isempty(tline), continue; end

    tline = strtrim(tline);
    if tline(1) == '!', continue; end
    
    break

end

end % GetNextLine

function S = sParmsBlank

    S = struct(...
        'Title', '' ...
        ,'strParms', '' ...
        ,'DataType', '' ...
        ,'ncols', 0 ...
        ,'nrows', 0 ...
        ,'DataValue', '' ...
        ,'wavelength', 0 ...
        ,'datascale_wave', 0 ...
        ,'nda', NaN ...
        );
    
end % sParmsBlank