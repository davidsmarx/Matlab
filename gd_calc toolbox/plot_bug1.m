function plot_bug
% This function illustrates an unresolved bug in MATLAB's 3-D plotting
% facility (as of Version 7.0.1, R14 SP2), which may affect the results of
% gdc_plot.m

for z=[.1 .05]
    figure
    axis([-1 1 -1 1 -1 1]);
    axis image
    xlabel('x');
    ylabel('y');
    zlabel('z');
    hold on;
    patch_([1;1;-1;-1],[-1;1;1;-1],[0;0;0;0],[1 0 0],'FaceAlpha',0.8);
    patch_([1;0;0;1],[-.5;.5;.5;-.5],[z;z;1;1],[0 1 0],'FaceAlpha',0.8);
    patch_([0;-1;-1;0],[-.5;.5;.5;-.5],[z;z;1;1],[0 1 0],'FaceAlpha',0.8);
end

function patch_(x,y,z,c,s,a);
if 1
    patch(x,y,z,c,s,a);
else
d=0.2;
for j=0:d:1-d/2
    for k=0:d:1-d/2
        patch(...
            [x(1)+j*(x(2)-x(1))+k*(x(4)-x(1));x(1)+(j+d)*(x(2)-x(1))+k*(x(4)-x(1));x(1)+(j+d)*(x(2)-x(1))+(k+d)*(x(4)-x(1));x(1)+j*(x(2)-x(1))+(k+d)*(x(4)-x(1))],...
            [y(1)+j*(y(2)-y(1))+k*(y(4)-y(1));y(1)+(j+d)*(y(2)-y(1))+k*(y(4)-y(1));y(1)+(j+d)*(y(2)-y(1))+(k+d)*(y(4)-y(1));y(1)+j*(y(2)-y(1))+(k+d)*(y(4)-y(1))],...
            [z(1)+j*(z(2)-z(1))+k*(z(4)-z(1));z(1)+(j+d)*(z(2)-z(1))+k*(z(4)-z(1));z(1)+(j+d)*(z(2)-z(1))+(k+d)*(z(4)-z(1));z(1)+j*(z(2)-z(1))+(k+d)*(z(4)-z(1))],...
            c,s,a,'EdgeColor','none');
    end
end
patch(x,y,z,c,'FaceColor','none');
end
