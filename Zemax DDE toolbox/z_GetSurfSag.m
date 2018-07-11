function [A, headerinfo] = z_GetSurfSag(zchan,textfilename)
% A = z_GetSurfSag(zchan,textfilename)
%
% zchan = dde interface
% textfilename (default = 'surfsag.txt') 

% ZEMAX performs the calculation with the latest saved settings for the
% requested analysis tool. Before calling this function, first check the
% current settings in ZEMAX, modify as necessary, and save.

% parse input arguments and set defaults
if nargin == 0, error('usage: A = z_GetSurfSag(zchan,textfilename)'); end

if ~exist('textfilename','var') || isempty(textfilename),
    textfilename = [pwd '\surfsag.txt']; 
end
%if exist(textfilename,'file'), delete(textfilename); end
typecode = 'Srs'; % default

% dummy file for zemax settings
settingsfilename = [pwd '\testsettings'];

% make the call
cmdstr = ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', 0'];
retval = ddereq(zchan,cmdstr,[1 1]);
%if ~strcmp(retval(1:2),'OK'), disp(['warning: dde return value: ', retval]); end

% wait until zemax is ready
while ~ddereq(zchan,'GetVersion',[1 1]), end

% read the data from the text file
[A, headerinfo] = ReadSurfSag(textfilename);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, headerinfo] = ReadSurfSag(textfilename)

unitsdefinitions;

fid = fopen(textfilename,'rt');
if fid == -1, error(['error opening text file ' textfilename]); end

% search for next data set
while ~feof(fid),
    ltmp = fgetl(fid);
    
    if ~isempty(ltmp) && ~isequal(ltmp, -1),
        wordparse = textscan(ltmp,'%s');
    
        switch lower(wordparse{1}{1}), % first word
            case 'file',
                atmp = textscan(ltmp,'File : %s','whitespace','\t');
                headerinfo.Lensfile = atmp{1}{1};
            case 'title:',

            case 'date',
                atmp = textscan(ltmp,'Date : %s','whitespace','\t');
                headerinfo.Date = atmp{1}{1};
                
            case 'configuration'
                atmp = textscan(ltmp,'Configuration %f');
                headerinfo.Configuration = atmp{1};
                
            case 'surface',
                atmp = textscan(ltmp,'Surface: %s','whitespace','\t');
                if ~isempty(atmp{1}),
                    headerinfo.Surface = atmp{1}{1};
                end
                
            case 'units'
                
            case 'center'
                atmp = textscan(ltmp,'Center point is: %f, %f');
                headerinfo.Centerpoint = [atmp{:}];
                
                % read the data
                ltmp = fgetl(fid); % blank line
                atmp = textscan(fid,'%f');
                A = reshape(atmp{1}, 2*headerinfo.Centerpoint - [1 1]);
                
            otherwise,
                % do nothing

        end
    end
end

