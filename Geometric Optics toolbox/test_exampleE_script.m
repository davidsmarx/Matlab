t = [15];
c = [0.02 -0.02];
nD = [1 1.5 1];
V = [inf 62.5 inf];

% V = (nD-1)/(nF-nC) = (nD-1)/Delta_n
dn = (nD - 1)./V;
n = [nD+dn/2; nD-dn/2]

[TSC, CC, TAchC,TchC,LAchC, EFL, li, hi] = paraxaberrations(t,c,n,20,0.1,0,0.1)