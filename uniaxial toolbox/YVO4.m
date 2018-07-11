lam = 1.55;
EFL = 1.9; % [mm] effective focal length of collimating lens
AOI = 2.3*pi/180; % desired angle of incidence for each ray at the collimating lens

no = sqrt( 3.77834 -0.0108133.*(lam.^2) + 0.069736./(lam.^2 -0.04724) );
ne = sqrt( 4.59905 -0.0122676.*(lam.^2) + 0.110534./(lam.^2 -0.04813) );

no = no + 8.5e-4;
ne = ne + 3.0e-4;

% wallaston prism:
th1 = [0:30]*pi/180; % angle of interface
% angle of rays w.r.t. interface normal:
thoep= asin( (ne./no).*sin(th1) );
theop= asin( (no./ne).*sin(th1) );
%fprintf('thoep = %f  theop = %f\n',thoep*180/pi,theop*180/pi);

% angle of rays w.r.t. optic axis:
thoe = th1 - thoep;
theo = th1 - theop;
%fprintf('thoe = %f  theo = %f\n',thoe*180/pi,theo*180/pi);

% exit angle from prism to air
thoea= asin( no.*sin(thoe) );
theoa= asin( ne.*sin(theo) );

figure, plot(th1*180/pi,[thoea(:) theoa(:)]*180/pi), legend('o-e','e-o'), grid
figure, plot(th1*180/pi,(thoea - theoa)*180/pi), grid

th11 = interp1(abs(thoea - theoa),th1,2*AOI); % find interface angle for desired AOI

fprintf('for AOI = %f deg, interface angle = %f deg\n',...
   AOI*180/pi,th11*180/pi);
fprintf('fiber separation = %f mm\n',EFL*tan(2*AOI));
