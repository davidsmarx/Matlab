T = 2048;
N = 20;

y = markov(T,N,1);

ym= mean(y,2);
ys= std(y,0,2);

h = repmat(hanning(T),1,N);
yf= fft(y.*h);
yp= sum(abs(yf).^2,2);
size(yp)

figure(1), plot(1:T,y)
figure(2), plot(1:T,ym)
figure(3), plot(1:T,ys)
figure(4), plot(1:T,abs(ym./ys))
figure(5), plot(1:T,yp)
