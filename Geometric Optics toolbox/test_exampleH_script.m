t = [8.74 11.05 2.78 7.63 9.54];
c = 1./[40.94 inf -55.65 39.75 107.56 -43.33];
nD = [1 1.617 1 1.649 1 1.617 1];
V = [inf 55.0 inf 33.8 inf 55.0 inf];
% V = (nD-1)/(nF-nC) = (nD-1)/Delta_n

dn = (nD - 1)./V;
n = [nD+dn/2; nD-dn/2]

[TSC, CC, TAchC,TchC,LAchC, EFL, li, hi] = paraxaberrations(t,c,n,18.5,0,-6.3,0.25)