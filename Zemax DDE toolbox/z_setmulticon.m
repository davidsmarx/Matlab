function retstr = z_setmulticon(zchan,nconf,row,value,solvetype,pickuprow,pickupconfig,scale,offset)
% status =
%    z_setmulticon(zchan,nconf,row,value,solvetype,pickuprow,pickupconfig,scale,offset)
% 
% set values in the multi-configuration editor
%
% zchan = dde handle
% nconf = configuration number
% row = row # in the multi-configuration editor
% value = new value
% solvetype = 0->fixed, 1->variable, 2->pickup, 3->thermal pickup
% if status = 2 or 3,
%   pickuprow and pickupconfig are source data for pickup
%   scale and offset apply to pickup source

ddestr = [
    'SetMulticon, '...
    num2str(nconf) ', '...
    num2str(row) ', '...
    num2str(value) ', '...
    num2str(solvetype) ', '...
    num2str(pickuprow) ', '...
    num2str(pickupconfig) ', '...
    num2str(scale) ', '...
    num2str(offset)
    ];

retstr = ddereq(zchan,ddestr,[1 1]);


return

