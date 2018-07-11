tic;
C = 2.9979e2;  % [um / 1e-12s]
f0 = 200; % [THz] center frequency
refl = C/f0; % [um] reference wavelength

[d, n] = get_thinfilmfilter(0,refl);

theta = 0.0;
freq = linspace(f0-0.2,f0+0.2);
lam  = C./freq;

[R, T] = filter_response(n,d,theta,lam);

phase = unwrap(angle(T));
phase = phase - mean(phase);
grvel = gradient(phase,freq(2)-freq(1));
chrom = gradient(grvel,freq(2)-freq(1));

figure, subplot(2,1,1), plot(freq-f0,decibel(T)), grid,...
   subplot(2,1,2), plot(freq-f0,phase/pi),grid
figure, plot(freq-f0,grvel), grid, xlabel('frequency offset [THz]'), ylabel('group velocity')
figure, plot(freq-f0,chrom), grid, ylabel('chromatic dispersion')

disp(['time: ' num2str(toc/60) ' minutes']);
