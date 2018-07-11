tic;
C = 2.99792458e2;
refl = C./191; % um reference wavelength

[d, n] = get_thinfilmfilter(0,refl);

%freq = C./refl + [0 1]; %[0:0.05:3.8];
lam  = refl; %C./freq;

theta = [0:19]*pi/180;

[Rte, Tte] = TFfilter_response(n,d,theta,lam,0); % size(R) = size(T) = [length(lam), length(theta)]
[Rtm, Ttm] = TFfilter_response(n,d,theta,lam,1); % size(R) = size(T) = [length(lam), length(theta)]

disp(toc/60);
return

figure, plot(lam,10*log10(Rplot),lam,10*log10(Tplot)),...
   grid, %set(gca,'ylim',[-20 0])

figure, plot(thetalist*180/pi,-1000*(refl-CWL),'b',thetalist*180/pi,-1000*(refl-CWL),'--k'), grid,...
   xlabel('Angle of Incidence [deg]'), ylabel('Wavelength blue shift [nm]'),...
   set(gca,'xlim',[1.8 2.9]), set(gca,'ylim',[-1.0 0.0]),...
   legend('200 GHz','100 GHz')

figure, plot(1000*(refl-CWL),thetalist*180/pi), grid,...
   ylabel('Angle of Incidence [deg]'), xlabel('Wavelength blue shift [nm]')

fprintf('%f minutes\n',toc/60);

return

% take gradient to find band edges and estimate center wavelength
dT = gradient(T(:));
[ml, il] = max(dT);
[fm, xl] = findpeak(dT(il-1),dT(il),dT(il+1));
[mr, ir] = min(dT);
[fm, xr] = findpeak(dT(ir-1),dT(ir),dT(ir+1));
cwl = 0.5*(lam(il)+xl*dlam + lam(ir)+xr*dlam);

% script for adding random layer thickness errors
Ntrials = 0;
R = zeros(length(lam),Ntrials);
T = zeros(length(lam),Ntrials);
for ii = 1:Ntrials,
   drand = d + (1.0e-4*0.25*refl./real(n(2:end-1))).*randn(size(d));
   for i = 1:length(lam),
      [r t rr tt] = thin_film_filter_2(n,drand,theta,lam(i),0);
      %[r t rr tt] = thin_film_filter(n(1),n(2:end),d,theta,lam(i),0);
      RR(i) = r(1);
      TT(i) = t(end);
   end
      
   R(:,ii) = RR(:);
   T(:,ii) = TT(:);
end

Rplotrand = [Rplotrand mean(R,2)];
Tplotrand = [Tplotrand mean(T,2)];

figure, plot(lam,10*log10(Tplot),'b',lam,10*log10(Tplotrand),'r'), grid, set(gca,'ylim',[-1 0]),...
   legend('pure','noisy'), xlabel('wavelength [\mum]'), ylabel('Response [dB]')
