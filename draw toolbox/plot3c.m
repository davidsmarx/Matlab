function hc = plot3c(x,y,z,v,marker,vlim,props)
%FUNCTION PLOT3C(X,Y,Z,V,'MARKER',vlim,props) plots the values in vector v colour coded
% at the positions specified by the vectors x, y, and z in a 3-D axis
% system. A colourbar is added on the right side of the figure.
%
% The colorbar strectches from the minimum value of v to its
% maximum in 9 steps (10 values).
%
% marker is optional to define the marker being used. The
% default is a point. To use a different marker (such as circles, ...) send
% its symbol to the function (which must be enclosed in '; see example).
%
% vlim is optional = [vmin vmax]
% hc = plot3c(...) returns the handle to the colorbar.
%
% props is optional = cell array of line property pairs, e.g.
%  props = {'MarkerSize',1}
%
% This function is an extension of PLOTC.
%
% Example:
% The seismic P-velocity (v) depends on the three parameters porosity (por) and the
% bulk moduli of the saturating fluid (kf) and the elastic frame (kd). To plot the
% velocity data as a function of these three parameters use (assuming that
% all data are given in vectors):
%
% plot3c(por,kd,kf,v,'d','Velocity')
%
% Uli Theune, University of Alberta, 2004
% modified by D. Marx, 12-04-2006

% initialize outputs
if nargout > 0, hc = []; end

%delete(gca)
if nargin < 3
    error('usage');
end
if nargin == 3,
    v = z;
end
if nargin <5
    marker='.';
end

if ~exist('props','var'),
    props = {};
end

map=colormap;
clen = length(map);
if exist('vlim','var') & ~isempty(vlim),
    miv = min(vlim);
    mav = max(vlim);
else,
    miv=min(v(:));
    mav=max(v(:));
end
eps = 1e-3*(mav-miv);
edges = (miv:(mav+eps-miv)/clen:mav+eps); % length(edges) should be clen+1
if isempty(edges),
    keyboard;
end

% Plot the points
i = 1;
in = v>=edges(i) & v<edges(i+1);
plot3(x(in),y(in),z(in),marker,'color',map(i,:),'markerfacecolor',map(i,:),props{:})
hold on
for i=2:length(edges)-1,
    in = v>=edges(i) & v<edges(i+1);
    hp = plot3(x(in),y(in),z(in),marker,'color',map(i,:),'markerfacecolor',map(i,:),props{:});
end
hold off

if mav <= miv
    miv = -1;
    mav = 1;
end
set(gca,'clim',[miv mav]);

% Re-format the colorbar
%hc=colorbar;
%hc  = [];

view(3)