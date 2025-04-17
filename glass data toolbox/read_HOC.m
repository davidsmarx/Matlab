function [A, unitstr] = read_HOC(fname)
% [A, unitstr] = read_HOC(fname)
% 
% fname = name of ascii file

%HOCDIR = 'C:\Tamar\Matlab\glass data toolbox\Handbook of Optical Constants\ascii tables\';
%HOCDIR = 'C:\Users\dmarx\OneDrive\Matlab\glass data toolbox\Handbook of Optical Constants\ascii tables\';
pn = regexp(userpath,'\','split');
HOCDIR = [pn{1} '\' pn{2} '\' pn{3} '\Documents\Matlab\glass data toolbox\Handbook of Optical Constants\ascii tables\'];

fid = fopen([HOCDIR fname],'rt');
if fid == -1, error(['cannot open file: ' HOCDIR fname]); end

lcnt = 1;
tline = fgetl(fid);
while ~strncmp(strtrim(tline),'eV',2),
   tline = fgetl(fid);
   lcnt = lcnt + 1;
end

unitstr = textscan(tline,'%s');
unitstr = unitstr{1};
Nunits = length(unitstr);

% next line is '-----------'
tline = fgetl(fid);
lcnt = lcnt + 2; 

A = textscan(fid,'%f','delimiter',':');

A = reshape(A{1},Nunits,length(A{1})/Nunits).';

fclose(fid);

% %A = textread([HOCDIR fname],'','headerlines',lcnt,'delimiter',':');
% A = dlmread([HOCDIR fname],':',lcnt,1);

end
