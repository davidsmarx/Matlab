function sMerged = UpdateStruct(s1, s2, varargin)
% sMerged = UpdateStruct(s1, s2, options)
% 
% update s1 with fields from s2
% if s1 and s2 have overlapping fieldnames, values from s2
%
% options:
%   'intersect' -- only update fields in s1, i.e. no new fields from s2
%
% taken from:
%   https://stackoverflow.com/questions/15245167/update-struct-via-another-struct-in-matlab

% overlapping field names
fnames_overlap = intersect(fieldnames(s1), fieldnames(s2));

if contains(varargin, 'intersect')
    % only update existing fields in s1
    sMerged = s1;
    for ii = 1:length(fnames_overlap)
        fntmp = fnames_overlap{ii};
        sMerged.(fntmp) = s2.(fntmp);
    end
    
    return
end


% remove overlapping fields from s1
sMerged = rmfield(s1, fnames_overlap);

% merge all remaining unique names
fnames = [fieldnames(sMerged); fieldnames(s2)];

% merge sMerged and s2
sMerged = cell2struct([struct2cell(sMerged); struct2cell(s2)], fnames, 1);

