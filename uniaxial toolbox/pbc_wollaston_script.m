n1 = 1.0;
%no = 1.54426; ne = 1.55335; % calcite
no = 1.94657144; ne = 2.15093137; % YVO4
aoi = 0.0; %2.3614*pi/180;
thw = 9.05*pi/180;
%aoi = 30*pi/180;
k1 = [0 sin(aoi) cos(aoi)]'; % for aoi in the x-y plane
%k1 = [0 sin(aoi) cos(aoi)]'; % for aoi in the y-z plane
%k1 = [0 sin(aoi) cos(aoi)]'; % copied from ZEMAX
%ei = [0 cos(aoi) sin(aoi)]'; % direction of incident E field
ei = [1 0 0]';

% first wedge:
c = [0 1 0]'; % crystal axis
n = [0 0 1]'; % surface normal
%
[ke, se, me, ko] = uniaxial_raytracing_iso2uni(k1, n, c, n1, no, ne);

% exit first wedge:
n = [0 sin(thw) cos(thw)]';
% ext-beam:
kea = snells_law(ke,n,me,1.0); % ray vector in air

% ord-beam:
koa = snells_law(ko,n,no,1.0);

% second wedge:
n = [0 sin(thw) cos(thw)]';
c = [1 0 0]';
% o-in
[koe, soe, moe, koo] = uniaxial_raytracing(ko,n,c,no,no,ne);

% e-in
[kee, see, mee, keo] = uniaxial_raytracing(ke,n,c,me,no,ne);
[eeo, eee, eer, teo, tee, re] = uniaxial_refraction(ke,keo,kee,see,n,[me;no;ne],c,ee);

% output of second wedge:
n = [0 0 1]';
c = [0 0 1]'; % shouldn't matter since output medium is air
%oo-in
[kooe, sooe, mooe, kooa] = uniaxial_raytracing(koo,n,c,no,n1,n1);
[eooa, eooe, eoor, tooa, tooe, roo] = uniaxial_refraction(koo,kooa,kooe,sooe,n,[no;n1;n1],c,eoo);
%oe-in
[koee, soee, moee, koea] = uniaxial_raytracing(koe,n,c,moe,n1,n1);
[eoea, eoee, eoer, toea, toee, roe] = uniaxial_refraction(koe,koea,koee,soee,n,[moe;n1;n1],c,eoe);
%eo-in
[keoe, seoe, meoe, keoa] = uniaxial_raytracing(keo,n,c,no,n1,n1);
[eeoa, eeoe, eeor, teoa, teoe, reo] = uniaxial_refraction(keo,keoa,keoe,seoe,n,[no;n1;n1],c,eeo);
%ee-in
[keee, seee, meee, keea] = uniaxial_raytracing(kee,n,c,mee,n1,n1);
[eeea, eeee, eeer, teea, teee, ree] = uniaxial_refraction(kee,keea,keee,seee,n,[mee;n1;n1],c,eee);

return

