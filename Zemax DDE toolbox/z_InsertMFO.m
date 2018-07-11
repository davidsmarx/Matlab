function retval = z_InsertMFO(zchan,noperand)
% retval = z_InsertObject(zchan,noperand)
%
% nsurf: surface
% noperand: line number of the operand to insert
% This item inserts a new optimization operand in the merit function editor. The syntax is:
% InsertMFO, operand
% The operand argument is an integer between 1 and the current number of operands plus 1, inclusive. The
% return value is the new number of operands. See also DeleteMFOobject number

cmdstr = ['InsertMFO,' num2str(noperand)];

retstr = ddereq(zchan,cmdstr,[1 1]);
retval = str2double(retstr);

return


