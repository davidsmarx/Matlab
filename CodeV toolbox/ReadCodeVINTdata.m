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
            sParms.ncols = str2double(wordparse{2});
            sParms.nrows = str2double(wordparse{3});
            wordparse = wordparse(4:end);
            
        case 'wfr'
            sParms.DataType = 'wavefront';
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
data = flipud(data.'); %
data(data==sParms.nda) = nan;
map = 2*pi*data/sParms.datascale_wave; % return map in radians

fclose(fid);

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