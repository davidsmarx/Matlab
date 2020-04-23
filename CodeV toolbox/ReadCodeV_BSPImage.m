function [Im, sParms] = ReadCodeV_BSPImage(fn, varargin)
% [Im, sParms] = ReadCodeV_BSPImage(fn, options)
%
% options
%   bDisplay = CheckOption('display', false, varargin{:});
%   titlestr = CheckOption('title', pwd2titlestr(fn), varargin{:});
%
% return
%   Im = complex field

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
    
    hfig = figure_mxn(1,2)
    hax(1) = subplot(1,2,1);
    imageschcit(x/U.MM, y/U.MM, abs(Im))
    xlabel('X (mm)'), ylabel('Y (mm)')
    title(titlestr, 'fontsize',14)
    colorbartitle('|Field|')

    hax(2) = subplot(1,2,2);
    imageschcit(x/U.MM, y/U.MM, angle(Im)/pi)
    xlabel('X (mm)'), ylabel('Y (mm)')
    colorbartitle('Field Phase (\pi rad)')
    
end % DisplayImage

function [Im, sParms] = ReadImTxt(fn)

U = CConstants;

fid = fopen(fn,'rt');
if isequal(fid,-1), error(['error opening file: ' fn]); end

% 
cnt = 0;

ltmp = fgetl(fid);

while ~feof(fid),
    if isempty(ltmp),
        ltmp = fgetl(fid);
        continue
    end
    
    wordparse = textscan(ltmp,'%s');
    
    switch wordparse{1}{1},
        case 'BSP'
            cnt = cnt + 1;
            
        case 'These'
            % These data represent the Intensity at surface 48
            % These data represent the Complex Field at surface 48
            sParms(cnt).DataType = wordparse{1}{5};
            sParms(cnt).AtSurface = wordparse{1}{end};
            
        case 'Zoom'
            A = textscan(ltmp,'Zoom Position: %d Field Number: %d Wavelength (nm.):	%f');
            [sParms(cnt).zoom, sParms(cnt).FieldNum, wvtmp] = deal(A{:});
            sParms(cnt).Wavelength = wvtmp * U.NM;
            
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
            sParms(cnt).dx = dx*U.MM;
            sParms(cnt).dy = dy*U.MM;
                        
        case 'Array'
            sParms(cnt).Nr = str2double(wordparse{1}{end-1});
            sParms(cnt).Nc = str2double(wordparse{1}{end});
            
            % Im data follows
            if strcmpi(sParms(cnt).DataType, 'Intensity'),
                Imtmp = fscanf(fid,'%f',[sParms(cnt).Nr sParms(cnt).Nc]);
            elseif strcmpi(sParms(cnt).DataType, 'Complex'),
                Imtmp = fscanf(fid,'%f',[2*sParms(cnt).Nc sParms(cnt).Nr]); % fscanf is in column order
                Imtmp = Imtmp(1:2:end,:) + 1i*Imtmp(2:2:end,:);
            end
            Im{cnt} = Imtmp.'; % back to row, column
            
        otherwise
            
    end % switch

    ltmp = fgetl(fid);
    
end % while ~feof

fclose(fid);

end % main

