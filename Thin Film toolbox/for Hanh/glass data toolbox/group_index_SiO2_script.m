% group index of amorphous SiO2 in visible spectrum
cls

% from HOC:
Table = [
    0.467816	1.46429
0.486133	1.46313
0.508582	1.46187
0.546074	1.46008
0.576959	1.45885
0.579065	1.45877
0.587561	1.45847
0.589262	1.45841
0.643847	1.45671
0.656272	1.45637
0.667815	1.45608
0.706519	1.45515
	
0.852111	1.45248
0.89435	1.45185
1.01398	1.45025
1.08297	1.44941
1.12866	1.44888
];
lam = Table(:,1)*UM;
n   = Table(:,2);

lami = linspace(500, 800)'*NM;

ng = index2groupindex(struct(...
    'wavelength', lam, 'index', n ...
    ), lami);

figure, plot(lam/NM, n, '-o', lami/NM, ng, '-'), grid on
xlabel('Wavelength (nm)')
ylabel('Index of Refraction')
legend('HOC','Group Index')
