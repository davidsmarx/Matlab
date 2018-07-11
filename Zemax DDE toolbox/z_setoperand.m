function opval = z_getoperand(zchan,row,column,opval)
% opval = z_setoperand(zchan,row,column,opval)

switch column
    case 1,
        cmdstr = sprintf('SetOperand,%d,%d,%s',row,column,opval);
        
    otherwise,
        cmdstr = sprintf('SetOperand,%d,%d,%f',row,column,opval);
end

opval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);


return
