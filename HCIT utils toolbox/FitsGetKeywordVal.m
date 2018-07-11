function [val, comment, ik] = FitsGetKeywordVal(Keywords, key)
% [val, comment, ik] = FitsGetKeywordVal(Keywords, key)
% 
% Keywords is cell array from finfo.PrimaryData (or Image) .Keywords
% key is string
%
% val is the keyword value
% comment is the comment (3rd column)
% ik is the row index where the keyword was found
%
% if key is not found, empty matrix is returned
% 
% if more than one key is found, val, comment are cell arrays

ik = find(strcmpi(key,Keywords(:,1)));

val = [];
if ~isempty(ik),
    if length(ik) == 1,
        val = Keywords{ik,2};
        comment = Keywords{ik,3};
        
        if val(end) == '&' && strcmp(Keywords{ik+1,1},'CONTINUE'),
            % continues
            stmp = strtrim(Keywords{ik+1,3});
            val = [val(1:end-1) stmp(2:end-2)];
        end
    else
        for ii = 1:length(ik)
            val{ii} = Keywords{ik(ii),2};
            comment{ii} = strtrim(Keywords{ik(ii),3});
        end
    end
end




