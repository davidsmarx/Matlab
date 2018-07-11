n1 = 1.0;
%no = 1.54426; ne = 1.55335; % calcite
no = 1.94657144; ne = 2.15093137; % YVO4
aoi = 0; %2.3614*pi/180;
k1 = [0 sin(aoi) cos(aoi)]'; % for aoi in the x-y plane
ei = [1 1 0]'; ei = ei./sqrt(ei'*ei);

% first wedge:
c = [sin(22.5*pi/180) cos(22.5*pi/180) 0]'; % crystal axis
n = [0 sin(10*pi/180) cos(10*pi/180)]'; % surface normal
%
[ke, se, me, ko] = uniaxial_raytracing_iso2uni(k1, n, c, n1, no, ne);
%[ke, se, me, ko, eetmp, eotmp] = uniaxial_raytracing_iso2uni_eig(k1, n, c, n1, no, ne);
[eo, ee, ers, erp, to, te, rs, rp] = uniaxial_match_iso2uni(k1,ko,ke,se,n,[n1;no;me],c,ei);
Eo = sqrt(no)*to*eo;
Ee = sqrt(ne)*te*ee;
Ers = sqrt(n1)*rs*ers;
Erp = sqrt(n1)*rp*erp;

% exit first wedge:
n = [0 sin(11.43*pi/180) cos(11.43*pi/180)]';
% ext-beam:
kea = snells_law(ke,n,me,1.0); % ray vector in air
[kere, sere, kero] = uniaxial_raytracing_uniref(ke,n,c,me,no,ne); % reflected e&o-rays

% ord-beam:
koa = snells_law(ko,n,no,1.0);
[kore, sore, koro] = uniaxial_raytracing_uniref(ko,n,c,no,no,ne);

