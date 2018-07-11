function [Keywords, ik] = FitsSetKeywordVal(Keywords, key, value)
% [Keywords, ik] = FitsSetKeywordVal(Keywords, key, value)
% 
% Keywords is cell array from finfo.PrimaryData (or Image) .Keywords
% key is string
%
% val is the keyword value
% ik is the row index where the keyword was found
%
% if key is not found, empty matrix is returned

ik = find(strcmpi(key,Keywords(:,1)));

if ~isempty(ik),
    Keywords{ik,2} = value;
end




