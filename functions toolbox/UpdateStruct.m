function sMerged = UpdateStruct(s1, s2)
% sMerged = UpdateStruct(s1, s2)
% 
% update s1 with fields from s2
% if s1 and s2 have overlapping fieldnames, values from s2
%
% taken from:
%   https://stackoverflow.com/questions/15245167/update-struct-via-another-struct-in-matlab

% remove overlapping fields from s1
sMerged = rmfield(s1, intersect(fieldnames(s1), fieldnames(s2)));

% merge all remaining unique names
fnames = [fieldnames(sMerged); fieldnames(s2)];

% merge sMerged and s2
sMerged = cell2struct([struct2cell(sMerged); struct2cell(s2)], fnames, 1);

