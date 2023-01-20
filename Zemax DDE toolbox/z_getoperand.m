function opval = z_getoperand(zchan,row,column)
% opval = z_getoperand(zchan,row,column)
% 
%   gets the operand from the Merit Function Editor

if nargin == 0,
    disp(['usage: opval = z_getoperand(zchan,row,column)']);
   return
end


cmdstr = sprintf('GetOperand,%d,%d',row,column);
opvalstr = ddereq(zchan,cmdstr,[1 1]);
    
if column == 1,
    % first column is operand type string
    % do nothing, just return the opval string
    opval = deblank(opvalstr);
    
else,
    % operand is a number
    opval = sscanf(opvalstr,'%f,',[1 inf]);
    
end


return
