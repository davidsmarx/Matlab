function [A, headerinfo] = z_GetEncEnergy(zchan,textfilename,type)
% A = z_GetEncEnergy(zchan,textfilename,type)
%
% zchan = dde interface
% textfilename (default = 'encenergy.txt') file name where
%         encircled energy results are stored.
% type = 'diff' (default) for diffraction encircled energy
%        'geom' for geometric encircled energy
%        'xse' for extended source encircled energy
%
% A = struct array with fields:
%   field = (x,y) coordinates of source field point  (Inf => diffraction
%      limit response)
%   refcoords = (x,y) coordinates of origin of encircled energy data
%   x = radius
%   y = encircled energy at radius x
%
% headerinfo = struct:
%   Lensfile
%   Date
%   Surface
%

% ZEMAX performs the calculation with the latest saved settings for the
% requested analysis tool. Before calling this function, first check the
% current settings in ZEMAX, modify as necessary, and save.

% parse input arguments and set defaults
if nargin == 0, error('usage: A = z_GetEncEnergy(zchan,textfilename)'); end

if ~exist('textfilename','var') | isempty(textfilename),
    textfilename = [pwd '\encenergy.txt']; 
end
if exist(textfilename,'file'), delete(textfilename); end
if exist('type','var'),
    switch lower(type),
        case {'diff','diffraction'},
            typecode = 'Enc';
        case {'geom','geometric'},
            typecode = 'Geo';
        case {'extended source','xse','extended'}
            typecode = 'Xse';
        otherwise,
            error('unrecognized type');
    end
else
    typecode = 'Enc'; % default
end

% dummy file for zemax settings
settingsfilename = [pwd '\testsettings'];

% make the call
cmdstr = ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', 0'];
retval = ddereq(zchan,cmdstr,[1 1]);
%if ~strcmp(retval(1:2),'OK'), disp(['warning: dde return value: ', retval]); end

% wait until zemax is ready
while ~ddereq(zchan,'GetVersion',[1 1]), end

% read the data from the text file
[A, headerinfo] = ReadEncEnergy(textfilename);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, headerinfo] = ReadEncEnergy(textfilename)

unitsdefinitions;

fid = fopen(textfilename,'rt');
if fid == -1, error(['error opening text file ' textfilename]); end

% n counts the number of data sets
n = 0;

% search for next data set
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
                
            case 'surface:',
                atmp = textscan(ltmp,'Surface: %s','whitespace','\t');
                headerinfo.Surface = atmp{1}{1};
                
            case 'diff',
                n = n+1;
                A(n).field = Inf;

                [A(n).refcoords, A(n).x, A(n).y] = readdata(fid,A(n));

            case 'field:'
                n = n+1;
                ftmp = textscan(ltmp,'Field: %f %s');
                switch ftmp{2}{1}, % units
                    case 'deg',
                        A(n).field = ftmp{1}*P;
                    case 'mm',
                        A(n).field = ftmp{1}*MM;
                    otherwise,
                        error(['unknown units: ' ftmp{2}{1}]);
                end

                [A(n).refcoords, A(n).x, A(n).y] = readdata(fid,A(n));

            otherwise,
                % do nothing

        end
    end
end

% catch,
%     fclose(fid);
%     error('ReadEncEnergy encountered an error');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [refcoords, x, y] = readdata(fid,A)

% convert units
unitsdefinitions;

dat = A;

tline = fgetl(fid); % empty line
tline = fgetl(fid); % reference coordinates
refcoords = sscanf(tline,'Reference Coordinates: %e %e')*MM;
tline = fgetl(fid); % column header

% read data until blank line
ii = 0;
tline = fgetl(fid);
while ~isempty(tline) & tline ~= -1,
    ii = ii + 1;
    a = sscanf(tline,'%f %f');
    x(ii) = a(1);
    y(ii) = a(2);
    tline = fgetl(fid);
end

x = x(:)*UM; % change to column vector
y = y(:);

return
