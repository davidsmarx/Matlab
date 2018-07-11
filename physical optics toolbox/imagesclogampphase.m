function haxout = imagesclogampphase(A,x,y,titlestr,varargin)
% hax = imagesclogampphase(A,x,y,titlestr)

if nargin < 2,
    [ny, nx] = size(A);
    x = (-nx/2:nx/2-1)';
    y = (-ny/2:ny/2-1)';
end

if ~exist('titlestr','var'), titlestr = ' '; end

pos = get(gcf,'position');
pos(3) = 1.75*pos(3);
set(gcf,'position',pos)

hax(1) = subplot(1,2,1);
% Aintensity = abs(A).^2;
% imagesc(x,y,abs(A)), axis image,...
%     title(pwd2titlestr(titlestr)), SetColorbar('Amp');
bAgt0 = abs(A) > 0;
minA = min(abs(A(bAgt0)));
imagesc(x,y,20*log10(abs(A)+1e-3*minA)), axis image,...
   title(pwd2titlestr(titlestr)), SetColorbar('dB');

hax(2) = subplot(1,2,2);
imagesc(x,y,angle(A)./pi), axis image, title(pwd2titlestr(titlestr)),...
    colorbartitle('\pi rad') %SetColorbar('\pi radians');

if nargout > 0,
    haxout = hax;
end

return

function htit = SetColorbar(titlestr)

h_cb = colorbar;
htit = get(h_cb,'title');
set(htit,'string',titlestr);

return