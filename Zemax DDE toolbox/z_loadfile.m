function status = z_loadfile(zchan, filename, surf_num)
% status = z_loadfile(zchan, filename, surf_num)
%
% filename = full path of zemax lens file
% surf_num (optional) append lens in file to current lens starting at
%    surface #surf_num
% status: 0 => success
%
% LoadFile also performs a GetUpdate

% check if filename has .zmx extension
if ~strcmpi(filename(end-3:end),'.zmx') && ~strcmpi(filename(end-3:end),'zmx"')
    filename = [filename '.zmx'];
end

% call update
cmdstr = ['LoadFile,' filename];
if exist('surf_num','var') & surf_num >= 0,
    cmdstr = [cmdstr ',' num2str(surf_num)];
end
retstr = ddereq(zchan,cmdstr,[1 1]);

status = str2num(retstr);

switch status,
    case 0,
        % no error
    case -999,
        % load file failed
        error(['Zemax Load File ' filename ' failed!']);
    case -1,
        % the update failed, cannot trace rays
        error(['Zemax Update Failed, Cannot Trace Rays']);
    otherwise
        error(['Unknown error: ' num2str(status)]);
end
        

return
