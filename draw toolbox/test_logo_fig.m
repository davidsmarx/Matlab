hf = tamarfigure;

data = load('vibesdat');
hs = surf(data.L1);

%%
%clf
data = load('vibesdat');
A = imread('tamarlogo.bmp');

hf = figure;
colordef(hf,'black');
hi = image(A);
ha = gca;
set(gca,'position',[0 0 0.2 0.2]);
axis off

%%

%shading interp
colormap(fliplr(pink));
ax = get(gca,'Position');
axes
h = surf(data.L1)
set(h,'EdgeColor',[0.3 0.3 0.6]);
axis off
