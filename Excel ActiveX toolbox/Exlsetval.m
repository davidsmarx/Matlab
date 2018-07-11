function rangehandle = Exlsetval(sheet,range,value,varargin)
% rangehandle = Exlsetval(sheet,range,value)
% rangehandle = Exlsetval(sheet,range,value,property,value)
% 
% range = cell array
% example:
%   Exlsetval(sheet, {[c1 r1],[c2 r2]}, values, 'NumberFormat', '0.00');
%   

if iscell(range),
	rangehandle = get(sheet,'range',range{:});
else
	rangehandle = get(sheet,'range',range);
end

set(rangehandle, 'value', value, varargin{:});

return
