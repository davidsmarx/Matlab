function [EFL, TAchC, TchC, TPC] =...
    ThinLensAberations(d,p,n,V,obj_z,obj_h,stop_h,stop_lens)
% [EFL, TAchC, TchC, TPC] = ThinLensAberations(d,p,stop)
%
% d = list of lens separations
% p = list of lens powers (p = 1/f)
% n = index of refraction of each lens. length(n) = length(p)
% V = Abbe V-number for each lens. length(V) = length(p)
% obj_z = distance from object to first lens
% obj_h = max object height
% stop_h = max stop height
% stop_lens = p(stop_lens) is the lens where the stop is located.
%
% length(d) must equal length(p). d(n) = distance from p(n) to p(n+1).
% d(end) = distance from last lens to image surface.
%
% EFL = effective focal length of system
% TAchC = transverse axial chromatic aberation
% TchC = lateral chromatic aberation
% TPC = transverse third order Petzval curvature (lateral displacement
% between where the ray hits the Petzval surface and where the same ray
% hits the image plane. Longitudinal Field Curvature = TPC/u(end). ZEMAX
% PETZ gives the radius of curvature of the Petzval surface, which is
% approximately yp(end)^2*u(end)/2/TPC.)

% trace a parallel ray to get focal length
[u, y] = thinlensraytrace(d,p,0,obj_h);
EFL = -y(1)/u(end);

% trace an axial ray
% for now put pupil at first lens, to do: modify to find axial ray when the
% pupil is not at the first lens.
[u, y] = thinlensraytrace(d,p,stop_h./obj_z,stop_h);

% determine image plane distance from d(end) and correct d(end)
image_z = y(end)./u(end);
d(end) = d(end) + image_z;

% trace a principal ray
% for now put pupil at first lens, to do: modify to find axial ray when the
% pupil is not at the first lens.
[up, yp] = thinlensraytrace(d,p,-obj_h./obj_z,0);

% Q factor for each lens (equation 10.8e)
Q = yp(1:end-1)./y(1:end-1);

% now calculate the aberations from the stop shift equations
tachc = -y(1:end-1).^2.*p./(V.*u(end));
tchc = -Q .* tachc;
tpc = 0.5*yp(end).^2.*p.*u(end)./n;

TAchC = sum(tachc);
TchC = sum(tchc);
TPC = sum(tpc);

EFL = simplify(EFL);
TAchC = simplify(TAchC);
TchC = simplify(TchC);
TPC = simplify(TPC);