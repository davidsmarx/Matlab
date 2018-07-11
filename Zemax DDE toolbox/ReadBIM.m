function bimimage = ReadBIM(filename)
% bimimage = ReadBIM(filename)
%
% filename = name of bim file to read (no path is required).
% BIM image files resulting from ZEMAX image analysis are automatically
% stored in the path C:\ZEMAX\IMAFiles

IMAdir = 'C:\ZEMAX\IMAFiles\';

if ~strcmp(filename(end-3:end),'.bim'), filename = [filename '.bim']; end

fid = fopen([IMAdir filename],'r');
nx = fread(fid,1,'int32');
ny = fread(fid,1,'int32');
bimimage = fread(fid,[ny nx],'double');

bimimage = rot90(bimimage);

fclose(fid);