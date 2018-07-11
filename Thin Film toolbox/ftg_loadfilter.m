function [val, Design] = ftg_loadfilter(filename)
% [val, Design] = ftg_loadfilter(filename)

Design = actxserver('Ftgdesign1.clsBasic');

val = invoke(Design,'FileOpen',filename);
if val ~= -1, error(['error opening ' filename]); end

if nargout < 2,
   release(Design);
end


return
