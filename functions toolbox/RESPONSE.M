function rcs = response(t,x_el,xmitter,object,freq,bw,el_width,velocity)
% rcs = response(t,x_el,xmitter,object,freq,bw,el_width,velocity)
% rcs: (1:Nt,1:Nl)
% x_el(1:Nl) = position of each receiver element
% xmitter(1:3) = (x,y,z) location of transmitter
% object(1:3) = (x,y,z) location of point object
% freq = Hz
% bw = Hz bandwidth
% el_width = [m]
% velocity = [m/s] of platform

if nargin ~= 8,
   error('usage: rcs = response(t,x_el,xmitter,object,freq,el_width,velocity)')
end

c  = 1500; % [m/s] speed of sound
Nt = length(t);
Nl = length(x_el);
rcs = zeros(Nt,Nl);
k = 2*pi*freq/c; 
fprintf('Nt = %d Nl = %d\n',Nt,Nl);

% distance and angle from transmitter to object:
Rto = sqrt( sum( (object - xmitter).^2 ) );
sto = (object(1) - xmitter(1))./Rto;  % = sin(theta)
fprintf('distance from transmitter to object = %f m\n',Rto);

% distance from object to each element, but first move the array
Ror = sqrt( (x_el+2*velocity*Rto/c - object(1)).^2 + object(2).^2 + object(3).^2 );
sor = (x_el+2*velocity*Rto/c - object(1)) ./ Ror;
fprintf('size(Ror) = %d x %d\n',size(Ror));

[R T] = meshgrid( Rto+Ror, t );
[Sor T] = meshgrid( sor, t );
fprintf('R: %d x %d\nT: %d x %d\nSor: %d x %d\n',...
   size(R),size(T),size(Sor));

% the pulse envelope
S = envelope(T-R/c,bw);
fprintf('size s = %d x %d\n',size(S));
rcs = sinc(k*el_width*sto/2).*sinc(k*el_width*Sor/2).*S.*exp(-j*k*R);
%rcs = S.*exp(j*k*R);

return

function s = envelope(t,bw)

c = 1500;
s = blackmanpuls(t,2/bw);
return

