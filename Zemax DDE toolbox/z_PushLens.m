function retstr = z_PushLens(zchan)
% retval = z_PushLens(zchan)

try
    retstr = ddereq(zchan,'PushLens, 1',[1 1],10000);

catch ME,
    warning('PushLens failed!');
    retstr = 'error';

end

return


