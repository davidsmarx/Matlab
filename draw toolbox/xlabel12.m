function hlabel = xlabel(varargin)
% xlabel('string')
% xlabel(fname)
% xlabel(...,'PropertyName',PropertyValue,...)
% xlabel(axes_handle,...)
% h = xlabel(...) 
%
% my default xlabel uses 'FontSize',12
%

if isempty(varargin{1}),
    varargin{1} = '';
end

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
htmp = xlabel(haxes,labelstring,'FontSize',12,proplist{:});
cd(curpath);

if nargout > 0, hlabel = htmp; end