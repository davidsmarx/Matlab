function a = range(b,dim)
% a = range(b)
% a = range(b,dim)
%
% a = max(b) - min(b)

switch nargin,
    case 1,
        % test if row or column array
        ss = size(b);
        if ss(1) == 1,
            dim = 2;
        else
            dim = 1;
        end
    case 2,
        % do nothing
    otherwise,
        error('usage: a = range(b,dim)');
end

a = max(b,[],dim) - min(b,[],dim);

