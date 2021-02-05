function [lam, val, sParms] = ReadOOspectrum(fn)
% [lam, val] = ReadOOspectrum(fn)

pn = './';
if nargin < 1 || isempty(fn),
    [fn, pn, fi] = uigetfile({'*.ProcSpec'});
end

%U = CConstants;
U.US = 1e-6; % microseconds
U.NM = 1e-9; % nanometers

fid = fopen([pn fn],'rt');

tline = fgetl(fid);
key = '>>>>>';
while ~feof(fid),
    
    if isempty(tline), tline = fgetl(fid); continue, end
    
    A = textscan(tline,'%s %s','delimiter',':');
    
    switch A{1}{1},
        case 'Number of Pixels in Processed Spectrum',
            sParms.Npix = str2double(A{2});
            
        case 'Integration Time (usec)',
            IntTime = sscanf(A{2}{1},'%f');
            sParms.IntTime = IntTime*U.US;
            
        case 'Correct for Electrical Dark',
            B = textscan(A{2}{1},'%s','delimiter',' ');
            if strcmp(B{1}{1},'No'),
                sParms.bDarkCorrection = false;
            else
                sParms.bDarkCorrection = true;
            end
                
            
        otherwise,
            if strcmp(tline(1:length(key)), key),
                % data follows
                A = textscan(fid,'%f %f',sParms.Npix,'delimiter','\t');
                [lam, val] = deal(A{:});
                lam = lam*U.NM;
                
                % next line should be '>>>>>End Processed Spectral
                % Data<<<<<'
    
                break;
            end
            
            
    end                        
    
    tline = fgetl(fid);
end

fclose(fid);
