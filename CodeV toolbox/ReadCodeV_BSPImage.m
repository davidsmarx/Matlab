function [Im, sParms] = ReadCodeV_BSPImage(fn, varargin)
% [Im, sParms] = ReadCodeV_BSPImage(fn, options)
%
% options
%   bDisplay = CheckOption('display', false, varargin{:});
%   titlestr = CheckOption('title', pwd2titlestr(fn), varargin{:});

% options
bDisplay = CheckOption('display', false, varargin{:});

if nargin == 0 || isempty(fn),
    [fn, pn] = uigetfile({'*.txt'});
    fn = [pn '/' fn];
end

[Im, sParms] = ReadImTxt(fn);

if bDisplay,
    [hfig, hax] = DisplayImage(Im, sParms, fn, varargin);
    sParms.hfig = hfig;
    sParms.hax = hax;
end

end % main

function [hfig, hax] = DisplayImage(Im, sParms, fn, varargin)

    U = CConstants;

    % options
    titlestr = CheckOption('title', pwd2titlestr(fn), varargin{:});
    
    [x, y] = CreateGrid(Im, sParms.dx);
    
    hfig = figure;
    imageschcit(x/U.MM, y/U.MM, Im)
    xlabel('X (mm)'), ylabel('Y (mm)')
    title(titlestr, 'fontsize',14)
    hax = gca;
    
end % DisplayImage

function [Im, sParms] = ReadImTxt(fn)

U = CConstants;

fid = fopen(fn,'rt');
if isequal(fid,-1), error(['error opening file: ' fn]); end

ltmp = fgetl(fid);

while ~feof(fid),
    if isempty(ltmp),
        ltmp = fgetl(fid);
        continue
    end
    
    wordparse = textscan(ltmp,'%s');
    
    switch wordparse{1}{1},
        
        case 'These'
            % These data represent the Intensity at surface 48
            % These data represent the Complex Field at surface 48
            sParms.DataType = wordparse{1}{5};
            sParms.AtSurface = wordparse{1}{end};
            
        case 'Zoom'
            A = textscan(ltmp,'Zoom Position: %d Field Number: %d Wavelength (nm.):	%f');
            [sParms.zoom, sParms.FieldNum, wvtmp] = deal(A{:});
            sParms.Wavelength = wvtmp * U.NM;
            
        case 'Defocus:'
            
        case 'Index'
            
        case 'Fraction'
            
        case 'Beam'
            
        case 'Coordinate'
            
        case 'Orientation'
            
        case 'Offset'
            
        case 'Grid'
            A = textscan(ltmp,'Grid spacing: %f %s %f');
            [dx, dy] = deal(A{[1 3]});
            units = A{2}{1}; % assum 'mm' for now
            sParms.dx = dx*U.MM;
            sParms.dy = dy*U.MM;
                        
        case 'Array'
            sParms.Nr = str2double(wordparse{1}{end-1});
            sParms.Nc = str2double(wordparse{1}{end});
            
            % Im data follows
            if strcmpi(sParms.DataType, 'Intensity'),
                Im = fscanf(fid,'%f',[sParms.Nr sParms.Nc]);
            elseif strcmpi(sParms.DataType, 'Complex'),
                Im = fscanf(fid,'%f',[2*sParms.Nc sParms.Nr]); % fscanf is in column order
                Im = Im(1:2:end,:) + 1i*Im(2:2:end,:);
            end
            Im = Im.'; % back to row, column
            
        otherwise
            
    end % switch

    ltmp = fgetl(fid);
    
end % while ~feof

fclose(fid);

end % main

