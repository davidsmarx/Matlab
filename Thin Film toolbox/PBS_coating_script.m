ns = 1.7;
nf = 1.6;
nl = 1.45;
theta = 45*pi/180;
lam = 1.55;
Nl = 7;

n = [ns repmat([nf nl],1,Nl) nf ns];
d = (0.25*lam./n(2:end-1)); % .* sqrt(1 - (sin(theta)*ns./n(2:end-1)).^2);

[Rs, T, rs, tt] = thin_film_filter_2(n,d,theta,lam,0);
[Rp, T, rp, tt] = thin_film_filter_2(n,d,theta,lam,1);

disp([Rs Rp]);
return

nflist = linspace(1.1,2.5);

for i = 1:length(nflist),
   nf = nflist(i);
   [R, T, rs(i), tt] = thin_film_filter_2([ns nf],[],theta,lam,0);
   [R, T, rp(i), tt] = thin_film_filter_2([ns nf],[],theta,lam,1);
end

figure, plot(nflist,rp,nflist,rs), grid

for i = 1:length(nflist),
   nf = nflist(i);
   [R, T, rs(i), tt] = thin_film_filter_2([nf ns],[],theta,lam,0);
   [R, T, rp(i), tt] = thin_film_filter_2([nf ns],[],theta,lam,1);
end

figure, plot(nflist,rp,nflist,rs), grid