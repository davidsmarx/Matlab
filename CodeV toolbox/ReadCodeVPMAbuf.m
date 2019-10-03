function [Sout] = ReadCodeVPMAbuf(fn, varargin)
% Sout = ReadCodeVIntensityTxt(fn)
%
% Sout(:) = struct(
%    pmap
%    ... parms ...
%    )
% for each field point and wavevelength

U = CConstants;

if nargin == 0 || isempty(fn),
    [fn, pn] = uigetfile({'*.txt'});
    fn = [pn '/' fn];
end

novalue = CheckOption('novalue', -99999, varargin{:});

iMap = 1; % count the number of maps

fid = fopen(strtrim(fn),'rt');
if isequal(fid,-1), error(['error opening file: ' fn]); end
Sout(iMap).Filename = fn;

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
            Sout(iMap).DataType = wordparse{1}{5};
            Sout(iMap).AtSurface = wordparse{1}{end};
            
        case 'Wavelength:'
            A = textscan(ltmp,'Wavelength: %f %s');
            [wvl, units] = deal(A{:});
            Sout(iMap).Wavelength = wvl * U.NM;
            
        case 'X'
            % field value or focual length
            
        case 'Y'
            % field value or focal length
            
        case 'Reference'
            A = textscan(ltmp,'Reference Sphere Radius: %f');
                        
        case 'Zoom'
            A = textscan(ltmp,'Zoom Position: %d Field Number: %d Wavelength Number: %d');
            [Sout(iMap).zoom, Sout(iMap).FieldNum, Sout(iMap).WaveNum] = deal(A{:});
            
        case 'Defocus:'
            
        case 'Array'
            A  = textscan(ltmp, 'Array Size: %d');
            N = A{1};
            % data follows
            A = textscan(fid, '%f', N*N);
            pmap = reshape(A{1}, [N N]);
            pmap(pmap == novalue) = NaN;
            Sout(iMap).pmap = pmap.'; % transpose to row, column
            
            
            % done with this map, increment
            iMap = iMap + 1;
        otherwise
            
    end % switch

    ltmp = fgetl(fid);
    
end % while ~feof

fclose(fid);

end % main

