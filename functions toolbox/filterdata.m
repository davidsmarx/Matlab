function varargout = filterdata(iuse,varargin)
% varargout = filterdata(iuse,varargin)
%
% [ap, bp, cp, ...] = filterdata(iuse, a, b, c, ...)
% [ap, bp, cp, ...] = filterdata(iuse, a, b, c, ..., holeval)
% ap = filterdata(iuse, a{:})
%
% iuse is a conditional, such as c > 0.5

switch length(varargin)
    
    case nargout,

        for ii = 1:nargout,
            varargout(ii) = {varargin{ii}(iuse)};
        end

    case nargout + 1,
    
        holeval = varargin{end};
        
        if isempty(holeval),
            for ii = 1:nargout,
                varargout(ii) = {varargin{ii}(iuse)};
            end

        else
            for ii = 1:nargout,
                vtmp = varargin{ii};
                vtmp(~iuse) = holeval;
                varargout(ii) = {vtmp};
            end
        end % if isempty(holeval)
        
    otherwise,
        switch nargout
            case 1,
                tmp = cell(size(varargin));
                for ii = 1:length(varargin),
                    tmp{ii} = varargin{ii}(iuse);
                end
                varargout(1) = {tmp};
            otherwise,
                error('usage');
        end
        
end


return
