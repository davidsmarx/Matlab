function [A, unitstr] = read_HOC(fname)
% [A, unitstr] = read_HOC(fname)
% 
% fname = name of ascii file

HOCDIR = 'C:\Tamar\Matlab\glass data toolbox\Handbook of Optical Constants\ascii tables\';

fid = fopen([HOCDIR fname],'rt');
if fid == -1, error(['cannot open file: ' HOCDIR fname]); end

lcnt = 1;
tline = fgetl(fid);
while isempty(tline) | ~strcmp(RemoveSpaces(tline),'eV'),
   tline = fgetl(fid);
   lcnt = lcnt + 1;
end
unitstr = tline;

lcnt = lcnt + 2; % next line is '-----------'

fclose(fid);

A = textread([HOCDIR fname],'','headerlines',lcnt,'delimiter',':');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bstr = RemoveSpaces(astr)

isp = isspace(astr);
bstr = astr(~isp);

bstr = bstr(1:2);

end