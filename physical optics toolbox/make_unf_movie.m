function M = make_unf_movie(unfname_inline,basename,list,plotfun,varargin)
% M = make_unf_move(unfname_inline,basename,list,plotfun,drawfun,args,...,drawfun,args,...)
%
% plotfun = [], @abs, @dba, or something else. The graph is plotfun(A),
% where A is the field amplitude from the unf files.

unitsdefinitions;

% validate inputs
if ~isa(plotfun,'function_handle'), plotfun = @abs; end

ifun = [];
if ~isempty(varargin),
    for ii = 1:length(varargin),
        if isa(varargin{ii},'function_handle'),
            ifun = [ifun ii];
        end
    end
end

hf = figure;

for ii = 1:length(list),
    
    filename = feval(unfname_inline,list(ii),basename);
    [ A, dx, dy, lambda, curv, x, y ] = ReadWavefrontUNF( filename );

    figure(hf), imagesc(x/MM,y/MM,feval(plotfun,A)), axis image
    
    % draw something
    for id = 1:length(ifun)-1,
%        drawfun = varargin{1};
        feval(varargin{ifun(id)},varargin{ifun(id)+1:ifun(id+1)-1});
    end
    feval(varargin{ifun(end)},varargin{ifun(end)+1:end});
    
    M(ii) = getframe;
    
end