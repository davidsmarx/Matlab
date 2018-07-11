function retval = z_DeleteMFO(zchan,noperand)
% retval = z_InsertObject(zchan,noperand)
%
% nsurf: surface
% noperand: line number of the operand to delete
%
% This item deletes an optimization operand in the merit function editor. The syntax is:
% DeleteMFO, operand
% The operand argument is an integer between 1 and the current number of operands, inclusive. The return
% value is the new number of operands. See also InsertMFO.

cmdstr = ['DeleteMFO,' num2str(noperand)];

retstr = ddereq(zchan,cmdstr,[1 1]);
retval = str2double(retstr);

return


