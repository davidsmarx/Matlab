clear all; close all;

% example of using thin_film_filter_2
%
% incident = air
% layer 1  = PMGI, index = 1.47, thickness = 500nm
% layer 2  = Ni, index from table, thickness = 50nm
% substrate = BK7, index = 1.5 for simplicity
% 
% angle of incidence = 10 degrees
%
% calculate response for wavelength from 522.5 to 577.5 nm

UM = 1e-6;
NM = 1e-9;
P  = pi/180;

% from Handbook of Optical Constants for Ni:
% wavelength (um), n, k
A = [
0.5166	1.71	3.06
0.5391	1.75	3.19
0.5636	1.8	3.33
0.5904	1.85	3.48
0.6199	1.93	3.65
0.6358	1.98	3.74
0.6526	2.02	3.82
0.6702	2.08	3.91
0.6888	2.14	4
0.7085	2.21	4.09
0.7293	2.28	4.18
0.7514	2.36	4.25
0.7749	2.43	4.31
0.7999	2.48	4.38
];
Ni_lam = A(1,:)*UM;
Ni_n   = A(2,:);
Ni_k   = A(3,:);

% calculate response
lamlist = linspace(522.5, 577.5)'*NM;
Nlam = length(lamlist);

% index of Ni for each wavelength:
nNi = interp1(Ni_lam, Ni_n, lamlist, 'linear');
kNi = interp1(Ni_lam, Ni_k, lamlist, 'linear');

[Rte, Tte, rrte, ttte, Rtm, Ttm, rrtm, tttm] = deal(zeros(Nlam,1));
for ilam = 1:Nlam,
    % inputs to thin_film_filter_2:
    % index [air PMGI Ni glass]
    n = [1.0 1.47 nNi(ilam)-1i*kNi(ilam) 1.5];
    % thickness [PMGI Ni]
    d = [500*NM 50*NM];
    % angle of incidence:
    theta = 10*P;
    
    [Rte(ilam), Tte(ilam), rrte(ilam), ttte(ilam)] = thin_film_filter_2(n,d,theta,lamlist(ilam),0);
    [Rtm(ilam), Ttm(ilam), rrtm(ilam), tttm(ilam)] = thin_film_filter_2(n,d,theta,lamlist(ilam),1);
    
end

figure,
subplot(2,1,1),
plot(lamlist/NM, [Tte Ttm]), grid
title('Transmitted Power')
xlabel('Wavelength (nm)'), ylabel('|t|^2')
legend('TE','TM')
subplot(2,1,2)
plot(lamlist/NM, [Rte Rtm]), grid
title('Reflected Power')
xlabel('Wavelength (nm)'), ylabel('|r|^2')

% plot transmitted amplitude and phase
figure, 
plotampphase(lamlist/NM, [ttte tttm], 'title', 'Transmission Coefficient', ...
    'xlabel','Wavelength (nm)', 'legend', 'TE','TM'),
