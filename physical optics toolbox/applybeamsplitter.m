function Ub = ApplyBeamsplitter(U0,rn,mainbs,direction)
% Ub = ApplyMainBeamsplitter(U0,rn,mainbs,direction)
%
% U0 = complex field vector
% rn = radius abscissae vector
% mainbs = struct describing the beamsplitter with fields:
%    D = vector of outer diameters defining regions within the beamsplitter
%    r,t = structs defining reflection and transmission properties. Fields
%    for r and t are:
%        c = reflection or transmission coefficient
%        apod = cell array, each element is a string, function handle, or
%            inline function, default = no apodization
%        p = cell array of cell arrays containing parameters relevant for
%            apod. feval(apod{1},rn,p{1}{:}) evaluates fun(rn,p1,p2,...)
% direction = 'reflection', or 'transmit'
%        

radius = 0.5*mainbs.D; % diameter to radius

switch lower(direction)
    case {'reflect', 'reflection'}
        Ub = BeamSplitter(U0,rn,radius,mainbs.r);
    case {'transmit', 'transmission'}
        Ub = BeamSplitter(U0,rn,radius,mainbs.t);
    otherwise
        error('unrecognized direction');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ub = BeamSplitter(U0,rn,rb,cc)
% U0 = complex field, rn = field abscissae, rb = radius of beamsplitter
% regions, cc = r or t struct

Ub = zeros(size(U0));

if length(rb) ~= length(cc), error('radii and coefficients not the same length'); end

cc = SetDefaults(cc);  % check valid struct properties

rb = [0; rb(:)];
for ir = 1:length(cc),
    nr = rn<rb(ir+1) & rn>=rb(ir);
    if ~isempty(cc(ir).apod), % apply apodization
        ccnr = GetApod(cc(ir),rn(nr));
    else,
        ccnr = cc(ir).c;
    end
        
    Ub(nr) = ccnr.*U0(nr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ccnr = GetApod(cc,rr)

ccnr = cc.c .* feval(cc.apod,rr,cc.p{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cc = SetDefaults(cc)

if ~isfield(cc,'apod'),
    [cc.apod] = deal([]); % default value is no apodization
    [cc.p]    = deal({});
end
