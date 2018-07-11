function retval = z_SaveLens(zchan, filename)
% retval = z_SaveLens(zchan, filename with full path)
% retval = z_SaveLens(zchan) calls uiputfile
%
% newly saved lens file is then updated, as with "GetUpdate"
%
% retval = 0 is success
% retval = -1, lens cannot be updated, ray tracing fails
% retval = -2 no file specified
% retval = -3 general error
% retval = -999, the file could not be saved

retval = -1;

if nargin == 1,
    [ftmp, pathname, fi] = uiputfile({'*.zmx'}, 'Save ZEMAX Lens File');
    if isequal(ftmp, 0),
        retval = -2;
        return;
    end
    filename = [pathname ftmp];
end

% DDE string limit:
if length(filename) > 240,
    error('filename is too long');
end

try
    cmdstr = ['SaveFile,' filename];
    retstr = ddereq(zchan, cmdstr, [1 1]);
    retval = str2double(retstr);

catch ME,
    warning('SaveFile failed!');
    retval = -3;

end

return


