function [value, rangehandle] = Exlgetval(sheet,range)
% [value, rangehandle] = Exlgetval(sheet,range)
% 
% range = cell array

if iscell(range),
	rangehandle = get(sheet,'range',range{:});
else
	rangehandle = get(sheet,'range',range);
end

ranval = rangehandle.value;

% convert cell array of cells to array of doubles if all elements are numbers
try,
   N = size(ranval);
   [cellval{1:prod(N)}] = deal(ranval{:});
   value = [cellval{:}];
   value = reshape(value,N);
catch,
   value = ranval;
end

return


%Range(Selection, Selection.End(xlDown)).Select