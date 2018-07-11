tic;
refl = 1;

%[d, n] = get_thinfilmfilter(0,refl);
%d = flipud(d); n = flipud(n);
n = 1 + rand(1,10); %[1 2.0 1.5]; % incident first
d = 0.25*refl./n(2:end-1);

freq = (1./refl)*linspace(0.999,1.001);
lam  = 1./freq;

for i = 1:length(lam),
   [R(i), T(i), rr(i), tt(i)] = thin_film_filter_2(n,d,0,lam(i),0);
end
toc

% now new one
[rb, ra, ta] = TFfilter2z(n,d,0,refl,0);
H = freqz(fliplr(rb),ra,pi*freq);
G = freqz(1,ta,pi*freq);
toc

%angle([H(:) rr(:)])/pi
%angle([G(:).*exp(-j*0.5*pi*freq(:)) tt(:)])/pi

figure, plot(freq,decibel([rr(:) H(:)])), grid
figure, plot(freq,decibel([tt(:) G(:)])), grid

figure, plot(freq,angle([rr(:) H(:)])/pi), grid
figure, plot(freq,angle([tt(:) G(:).*exp(-j*length(d)*0.5*pi*freq(:))])/pi), grid