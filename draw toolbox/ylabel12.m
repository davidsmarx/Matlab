function hlabel = ylabel(varargin)
% ylabel('string')
% ylabel(fname)
% ylabel(...,'PropertyName',PropertyValue,...)
% ylabel(axes_handle,...)
% h = ylabel(...) 
%
% my default ylabel uses 'FontSize',12
%

switch class(varargin{1}),
    case {'char', 'function_handle'}
        haxes = gca;
        labelstring = varargin{1};
        proplist = {varargin{2:end}};
    case 'double' % must be a handle
        haxes = varargin{1};
        labelstring = varargin{2};
        proplist = {varargin{3:end}};
    otherwise,
        error('usage');
end
     
curpath = pwd;
cd([matlabroot '\toolbox\matlab\graph2d']);
htmp = ylabel(haxes,labelstring,'FontSize',12,proplist{:});
cd(curpath);

if nargout > 0, hlabel = htmp; end