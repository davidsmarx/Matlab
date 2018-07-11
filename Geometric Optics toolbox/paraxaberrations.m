function [TSC,CC,TAchC,TchC,LAchC,EFL,li,hi] = paraxaberrations(varargin)
% [TSC,CC,TAchC,TchC,LAchC,EFL,li,hi] = paraxaberrations(t,c,n,obj_h,obj_d,pupil_h)
% [TSC,CC,TAchC,TchC,LAchC,EFL,li,hi] = paraxaberrations(t,c,n,ya,ua,yp,up)
% obj_h = object height
% obj_d = distance from object to first surface
% pup_h = pupil height
% ya, ua = axial ray height and slope to first surface
% yp, up = principle ray height and slope to first surface
% t  = list of thicknesses, length(d) = # surfaces - 1
% c  = list of surface curvatures, length(c) = # surfaces
% n  = 2 x length(c)+1 matrix, where n(1,:) are the index at the short
%      wavelength, and n(2,:) are at the long wavelength. If n is a vector,
%      chromatic aberrations are not calculated.
% 
% TSC = paraxial spherical aberration contribution from each surface
% CC  = paraxial coma contribution from each surface
% EFL = effective focal length = axial ray input height / axial ray image slope
% li  = distance from last surface to plane where axial ray crosses axis
% hi  = height of paraxial ray at image plane

[t,c,n,dn,ya,ua,yp,up] = ValidateInput(varargin{:});

% calculate EFL
af = [0; ya];
[uf,yf,iif] = paraxraytrace_surf(t,c,n,af);
EFL = ya./uf(end);

% axial ray
%u0a = p1./l1; y0a = p1;
aa = [ua; ya];
[ua,ya,iia] = paraxraytrace_surf(t,c,n,aa);
li = -ya(end)/ua(end);

% principle ray
%u0p = h0./l1; y0p = 0;
ap = [up; yp];
[up,yp,iip] = paraxraytrace_surf(t,c,n,ap);
hi = yp(end) + li*up(end);

invar = ap(2)*n(1)*aa(1) - aa(2)*n(1)*ap(1);
hcheck = invar./(n(end)*ua(end));

nn = n(1:end-1); np = n(2:end);

Ba = nn.*(np-nn).*ya.*(ua + iia)./(2.*np.*invar);
Bp = nn.*(np-nn).*yp.*(up + iip)./(2.*np.*invar);

TSC = Ba.*(iia.^2).*hi;
CC = Ba.*iia.*iip.*hi;

% do chromatic aberrations if dispersion is given
if ~isempty(dn),
    dnn = dn(1:end-1); dnp = dn(2:end);
    
    TAchC = -ya.*iia.*(dnn - nn.*dnp./np)./(np(end).*ua(end));
    TchC = -ya.*iip.*(dnn - nn.*dnp./np)./(np(end).*ua(end));
    LAchC = -TAchC./ua(end);
else,
    TAchC = [];
    TchC = [];
    LAchC = [];
end
    

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,c,n,dn,ya,ua,yp,up] = ValidateInput(varargin)
% [TSC,CC,EFL,li,hi] = paraxaberrations(t,c,n,obj_h,obj_d,pupil_h)
% [TSC,CC,EFL,li,hi] = paraxaberrations(t,c,n,ya,ua,yp,up)

switch length(varargin),
    case 6,
        obj_h = varargin{4};
        obj_d = varargin{5};
        pup_h = varargin{6};
        yp = 0;
        up = obj_h/obj_d;
        ya = pup_h;
        ua = pup_h/obj_d;
    case 7,
        ya = varargin{4};
        ua = varargin{5};
        yp = varargin{6};
        up = varargin{7};
    otherwise,
        error('usage: [TSC,CC,EFL,li,hi] = paraxaberrations(t,c,n,ya,ua,yp,up)');
end
t = varargin{1};
c = varargin{2};
n = varargin{3};

if length(t) ~= length(c)-1, error('paraxaberrations: invalid vector t'); end

% test n, if n is a matrix, then determine dispersion for chromatic
% aberrations
[Mdisp, Mmedia] = size(n);

if Mmedia ~= length(c)+1, error('paraxaberrations: invalid vector n'); end

switch Mdisp,
    case 1,
        dn = [];
    case 2,
        dn = n(1,:) - n(2,:);
        n = mean(n); % row vector of length Mmedia
    otherwise,
        error('invalid index argument');
end

return
   