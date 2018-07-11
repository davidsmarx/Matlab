function varargout = lsdir(strdirname)

if nargin > 0,
    d = dir(strdirname);
else,
    d = dir;
end

strdirlist = d([d.isdir]);

switch nargout,
    case 0,
        disp(strvcat(strdirlist.name));
    case 1,
        varargout = {strvcat(strdirlist.name)};
    otherwise,
        error('too many output arguments');
end

return

