% 
function [lam, val] = ReadOOspectrum(fn)
% [lam, val] = ReadOOspectrum(fn)

pn = './';
if nargin < 1 || isempty(fn),
    [fn, pn, fi] = uigetfile({'*.ProcSpec'});
end

fid = fopen([pn fn],'rt');

tline = fgetl(fid);
key = '>>>>>';
while ~feof(fid),
    
    if isempty(tline), tline = fgetl(fid); continue, end
    
    A = textscan(tline,'%s %s','delimiter',':');
    
    switch A{1}{1},
        case 'Number of Pixels in Processed Spectrum',
            Npix = str2double(A{2});
            
        otherwise,
            if strcmp(tline(1:length(key)), key),
                % data follows
                A = textscan(fid,'%f %f',Npix,'delimiter','\t');
                [lam, val] = deal(A{:});
                
                % next line should be '>>>>>End Processed Spectral
                % Data<<<<<'
    
                break;
            end
            
            
    end                        
    
    tline = fgetl(fid);
end

fclose(fid);
