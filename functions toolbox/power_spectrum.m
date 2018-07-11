function ps = power_spectrum(s,N)
% ps = power_spectrum(s,N)
% power spectrum estimation
% s = input sampled data
% N = number of frequency bins for estimated spectrum

[nr nc] = size(s);
if nc ~= 1,
   error('input must be a column vector')
end
if nr < 2*N,
   error('input must be a longer vector')
end

PS = zeros([N 1]);
%s = s - mean(s);
%normalize:
sumw = sum(hanning(N).^2);
den = 0;
for i = 1 : N/2 : nr-N,
   ss = s(i:i+N-1);
   ss = ss - mean(ss);
   sp = hanning(N).*ss;
   PS = PS + (abs(fft(sp)).^2)/N;
   den = den + sumw;
end

%PS(1) = PS(N);
ps = fftshift(PS);
ps = PS/den;