function [val, comment, ik] = FitsGetKeywordVal(Keywords, key)
% [val, comment, ik] = FitsGetKeywordVal(Keywords, key)
% 
% Keywords is cell array from finfo.PrimaryData (or Image) .Keywords
% key is string or a cell array of strings
%
% val is the keyword value or a cell array of values
% comment is the comment (3rd column) or a cell array of comments
% ik is the row index (or array of index) where the keyword was found
%
% if key is not found, empty matrix is returned
% 
% if more than one key is found, val, comment are cell arrays

switch class(key),
    case 'char'
        [val, comment, ik] = GetOneKeyVal(key);

    case 'cell'
        listKey = key;

        val = cell(length(listKey),1);
        comment = cell(length(listKey),1);
        ik = zeros(length(listKey),1);
        for ikey = 1:length(listKey),
            [vtmp, ctmp, itmp] = GetOneKeyVal(listKey{ikey});
            if ~isempty(vtmp),
                [val{ikey}, comment{ikey}, ik(ikey)] = deal(vtmp, ctmp, itmp);
            else
                [val{ikey}, comment{ikey}, ik(ikey)] = deal(NaN, ' ', -1);
            end
        end

    otherwise,
        error('usage: key must be char or cell array');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val, comment, ik] = GetOneKeyVal(key)
        ik = find(strcmpi(key,Keywords(:,1)));
        
        val = []; comment = [];
        if ~isempty(ik),
            if length(ik) == 1,
                val = Keywords{ik,2};
                comment = Keywords{ik,3};
                                
                if ~isempty(val) && val(end) == '&' && strcmp(Keywords{ik+1,1},'CONTINUE'),
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
        
    end % GetOneKeyVal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % main

