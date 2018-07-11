function [hf, ha] = tamarfigure(varargin)
% [hf, ha] = tamarfigure(propertyedefs)
%
% creates a new figure window with hf = figure(propertydefs{:})
% adds the Tamar logo and sets the figure properties to
% Tamar defaults
%
% ha = handle to a new axes

%clf
A = imread('tamarlogo.bmp');

figure(varargin{:});
colordef(gcf,'black');
hi = image(A);
xw = .7*.2;
yw = .7*.1;
set(gca,'position',[1-xw 1-yw xw yw]);
axis equal
axis off

% create a new axes for plotting
axes;

% output arguments
if nargout > 0,
    hf = gcf;
    ha = gca;
end

% %%
% 
% %shading interp
% colormap(fliplr(pink));
% ax = get(gca,'Position');
% axes
% h = surf(data.L1)
% set(h,'EdgeColor',[0.3 0.3 0.6]);
% axis off
