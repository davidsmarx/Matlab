function h_plot=gdc_plot(grating,param_index,pmt_display,x_limit,new_fig)
% gdc_plot
%
% Plot a biperiodic grating structure (for use with gdc).
%
% Syntax:
%
%   h_plot=gdc_plot(grating,param_index,pmt_display,x_limit);
%   Generate the plot in a new figure window.
%
%   h_plot=gdc_plot(grating,param_index,pmt_display,x_limit,false);
%   Generate the plot in an existing figure window.
%
% Documentation reference:
%   Grating Diffraction Calculator (GD-Calc TM)
%   Coupled-Wave Theory for Biperiodic DiffractionGratings
%   (GD-Calc.pdf)
%
% Notes:
%
%   If there is a very large number of plot elements gdc_plot will give the
%   user the opportunity to cancel the plot.
%
%   MATLAB's 3-D plotting facility is currently very buggy (as of Version
%   7.0.1, R14 SP2). Overlapping surface patches may not display properly
%   (especially if transparency is used), and line edges may show through
%   opaque surfaces. (plot_bug1.m demonstrates the bug.) The gdc_plot code
%   mitigates these problems somewhat by displaying a small clearance
%   (x_eps) between the grating's structural blocks, but until MathWorks
%   fixes the bug there may still be problems. When it gets fixed x_eps can
%   be set to a very small or zero value.
%
% inputs:
%
%   grating: same as first gdc input
%
%   param_index (size-[1,?] integer): multidimensional parameter indices
%   for selecting which parameter combination to display. param_index must
%   be compatible with parameter sizes, allowing for repmat extension of
%   singleton dimensions. (Run param_size=gdc(grating) to get parameter
%   sizes and ensure that 1<=param_index(j)<=param_size(j), except the
%   second relation need not hold if param_size(j)==1.)
%
%   pmt_display (size-[1,?] struct): display properties for optical
%   materials (permittivities). pmt_display must be size-matched to
%   grating.pmt, and comprises the following fields (all required):
%
%     pmt_display(m).name (string): legend name for material grating.pmt{m}
%     If the string is empty ('') the material will not be shown in the
%     legend.
%
%     pmt_display(m).color (size-[1,3] real, or []): RGB color for
%     displaying material grating.pmt{m}, or [] to suppress displaying the
%     material. Values must be in the range 0 <= RGB <= 1.
%
%     pmt_display(m).alpha (real): transparency value for material
%     grating.pmt{m} The value must be in the range 0 <= alpha <= 1. (0
%     means transparent; 1 means opaque.)
%
%   x_limit (size-[2,3] real): limits of display rectangle: x_limit =
%   [x1_min, x2_min, x3_min; x1_max, x2_max, x3_max], with x_limit(1,:) <
%   x_limit(2,:). (The grating substrate is the half-space x1<=0.)
%
%   new_fig (true or false, optional input, default=true): Set new_fig to
%   false to plot into the current figure; otherwise a new figure will be
%   created.
%
% output:
%
%   h_plot (handle or []): handle to plot window, or [] if user canceled
%   plot
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

h_plot=[];

% Validate inputs
if ~any(nargin==[4,5])
    error('Wrong number of input arguments.');
end
param_size=gdc(grating);
if ~isnumeric(param_index) || ~isreal(param_index) || ...
        ~isequal(param_index,fix(param_index)) || ...
        ndims(param_index)>2 || size(param_index,1)~=1
    error('Wrong data type or size (param_index).');
end
param_size(1,end+1:length(param_index))=1;
param_index(1,end+1:length(param_size))=1;
param_index(param_size==1 & param_index>1)=1;
if any(param_index<=0) || any(param_index>param_size)
    error('param_index is out of range.');
end
if ~isstruct(pmt_display) || ~isequal(size(pmt_display),size(grating.pmt))
    error('Wrong data type or size (pmt_display).');
end
if ~isfield(pmt_display,'name')
    error('Missing data field (pmt_display(m).name).');
end
if ~isfield(pmt_display,'color')
    error('Missing data field (pmt_display(m).color).');
end
if ~isfield(pmt_display,'alpha')
    error('Missing data field (pmt_display(m).alpha).');
end
for m=1:length(pmt_display)
    if ~ischar(pmt_display(m).name)
        error('Wrong data type (pmt_display(%d).name).',m);
    end
    if ~isequal(pmt_display(m).color,[])
        if ~isnumeric(pmt_display(m).color) || ...
                ~isreal(pmt_display(m).color) || ...
                ~isequal(size(pmt_display(m).color),[1,3])
            error('Wrong data type or size (pmt_display(%d).color).',m);
        end
        if min(pmt_display(m).color)<0 || max(pmt_display(m).color)>1
            error('pmt_display(%d).color must be in the range 0 to 1.',m);
        end
    end
    if ~isnumeric(pmt_display(m).alpha) || ...
            ~isreal(pmt_display(m).alpha) || ...
            length(pmt_display(m).alpha)~=1
        error('Wrong data type or size (pmt_display(%d).alpha).',m);
    end
    if pmt_display(m).alpha<0 || pmt_display(m).alpha>1
        error('pmt_display(%d).alpha must be in the range 0 to 1.',m);
    end
end
if ~isnumeric(x_limit) || ~isreal(x_limit) || ~isequal(size(x_limit),[2,3])
    error('Wrong data type or size (x_limit)');
end
if ~all(x_limit(1,:)<x_limit(2,:))
    error('x_limit must satisfy x_limit(1,:)<x_limit(2,:).')
end
if nargin<5
    new_fig=true;
end
if ~isequal(new_fig,true) && ~isequal(new_fig,false)
    error('new_fig must be true or false.');
end

% Warn user if there are many structure elements in the figure.
N=0;
for l1=1:length(grating.stratum)
    N=N+num_elements(grating.stratum{l1});
end
if N>100 && ~isequal(questdlg(['WARNING: There are ' num2str(N) ...
        ' structure elements. OK to proceed with plot?']),'Yes')
    return
end

% Initialize figure.
if new_fig
    h_plot=figure;
else
    h_plot=gcf;
    cla;
end
axis([...
    x_limit(1,2),x_limit(2,2),...
    x_limit(1,3),x_limit(2,3),...
    x_limit(1,1),x_limit(2,1)...
    ]);
axis image
xlabel('x_2');
ylabel('x_3');
zlabel('x_1');
hold on;

% Allocate patch handles for generating legend.
h=zeros(size(grating.pmt));

% x_eps = display clearance between structural blocks.
% (Note: The clearance is required because MATLAB's patch function does not
% reliably display intersecting or close parallel surface patches.)
x_eps=max(x_limit(2,:)-x_limit(1,:))*1.0e-3;

% x_max = lateral radius enclosing x_limit
x_max=sqrt(max(x_limit(:,2).^2)+max(x_limit(:,3).^2))+x_eps;

% Display the substrate.
m=grating.pmt_sub_index;
if ~isempty(pmt_display(m).color)
    x0=x_limit(1,:)+x_eps;
    edge1=[-x0(1),0,0];
    if edge1(1)>0
        edge2=[0,x_limit(2,2)-x_eps-x0(2),0];
        edge3=[0,0,x_limit(2,3)-x_eps-x0(3)];
        h_=fill_box(x0,edge1,edge2,edge3,x_limit,...
            pmt_display(m).color,pmt_display(m).alpha);
        if h_~=0
            h(m)=h_;
        end
    end
end

% Grating periods:
d21=param(grating.d21,param_index);
d31=param(grating.d31,param_index);
d22=param(grating.d22,param_index);
d32=param(grating.d32,param_index);

% Loop through strata.
x1_top=0;
dx2=0;
dx3=0;
for l1=1:length(grating.stratum)
    [x1_top,dx2,dx3,h]=plot_stratum(...
        x1_top,dx2,dx3,h,grating.stratum{l1},...
        param_index,pmt_display,x_limit,x_eps,x_max,d21,d31,d22,d32);
end

% Construct legend.
h_=[];
name={};
for m=1:length(h)
    if h(m)~=0 && ~isempty(pmt_display(m).name)
        h_(end+1)=h(m);
        name{end+1}=pmt_display(m).name;
    end
end
if ~isempty(h_)
    legend(h_,name{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function N=num_elements(stratum)
if stratum.type==0
    N=1;
elseif stratum.type==1
    N=length(stratum.stripe);
elseif stratum.type==2
    N=0;
    for l2=1:length(stratum.stripe)
        if stratum.stripe{l2}.type==0
            N=N+1;
        else
            N=N+length(stratum.stripe{l2}.block);
        end
    end
elseif stratum.type==3
    N=0;
elseif stratum.type==4
    N=0;
    for l1=1:length(stratum.stratum)
        N=N+num_elements(stratum.stratum{l1});
    end
    N=stratum.rep_count*N;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [x1_top,dx2,dx3,h]=plot_stratum(x1_top,dx2,dx3,h,stratum,...
    param_index,pmt_display,x_limit,x_eps,x_max,d21,d31,d22,d32)
h_fig=gcf;
drawnow;
% The next line is a bug work-around to make the patch graphics go to the
% correct window after the user responds Yes to the questdlg prompt.
% (The bug has been fixed as of MATLAB Version 7.1.0 (R14), SP3.) 
figure(h_fig);
if stratum.type==0
    % homogeneous stratum
    % (Only the top and bottom stratum surfaces are displayed; the interior
    % appears hollow.)
    x1_bot=x1_top;
    x1_top=x1_bot+param(stratum.thick,param_index);
    m=stratum.pmt_index;
    if ~isempty(pmt_display(m).color)
        x0=[x1_bot,x_limit(1,2)-x_eps,x_limit(1,3)-x_eps];
        edge1=[x1_top-x1_bot,0,0];
        if edge1(1)>2*x_eps
            x0(1)=x0(1)+x_eps;
            edge1(1)=edge1(1)-2*x_eps;
        else
            x0(1)=x0(1)+edge1(1)/4;
            edge1(1)=edge1(1)/2;
        end
        edge2=[0,x_limit(2,2)-x_limit(1,2)+2*x_eps,0];
        edge3=[0,0,x_limit(2,3)-x_limit(1,3)+2*x_eps];
        h_=fill_box(x0,edge1,edge2,edge3,x_limit,...
            pmt_display(m).color,pmt_display(m).alpha);
        if h_~=0
            h(m)=h_;
        end
    end
elseif any(stratum.type==[1,2])
    % periodic stratum
    % (Homogeneous stripes are displayed as hollow tubes; structural blocks
    % are displayed as hollow boxes.)
    x1_bot=x1_top;
    x1_top=x1_bot+param(stratum.thick,param_index);
    if stratum.type==1
        % uniperiodic stratum
        % Compute stratum period (GD-Calc.pdf, equations 3.22, 3.23).
        fs1=[stratum.h11,stratum.h12]/[d21,d22;d31,d32];
        [fs21,fs31]=deal(fs1(1),fs1(2));
        den=fs21^2+fs31^2;
        [ds21,ds31]=deal(fs21/den,fs31/den);
        % Define unit basis vector [e21,e31] perpendicular to the stratum's
        % stripes, and unit basis vector [e22,e32] parallel to the stripes.
        den=sqrt(ds21^2+ds31^2);
        e21=ds21/den;
        e31=ds31/den;
        e22=e31;
        e32=-e21;
        clear fs1 fs21 fs31 den
    else
        % biperiodic stratum
        % Compute stratum period (GD-Calc.pdf, equation 3.18).
        ds=[d21,d22;d31,d32]/...
            [stratum.h11,stratum.h12;stratum.h21,stratum.h22];
        [ds21,ds22,ds31,ds32]=deal(ds(1,1),ds(1,2),ds(2,1),ds(2,2));
        % Define unit basis vector [e21,e31] perpendicular to the stratum's
        % stripes, and unit basis vector [e22,e32] parallel to the stripes.
        den=sqrt(ds22^2+ds32^2);
        e22=ds22/den;
        e32=ds32/den;
        tmp=[ds21,ds31]-[ds21,ds31]*[e22;e32]*[e22,e32];
        tmp=tmp/sqrt(tmp(1)^2+tmp(2)^2);
        e21=tmp(1);
        e31=tmp(2);
        clear ds den tmp
    end
    % Loop through stripes.
    c1=param(stratum.stripe{end}.c1,param_index)-1;
    for l2=1:length(stratum.stripe)
        stripe=stratum.stripe{l2};
        c1_=c1;
        c1=param(stripe.c1,param_index);
        if stratum.type==1 || stripe.type==0
            % homogeneous stripe
            m=stripe.pmt_index;
            if ~isempty(pmt_display(m).color)
                x0=[x1_bot,c1*[ds21,ds31]+[dx2,dx3]];
                x0(2:3)=x0(2:3)*[e21;e31]*[e21,e31]-x_max*[e22,e32];
                edge1=[x1_top-x1_bot,0,0];
                edge2=[0,e22,e32]*2*x_max;
                edge3=[0,(c1_-c1)*[ds21,ds31]*[e21;e31]*[e21,e31]];
                if edge1(1)>2*x_eps
                    x0(1)=x0(1)+x_eps;
                    edge1(1)=edge1(1)-2*x_eps;
                else
                    x0(1)=x0(1)+edge1(1)/4;
                    edge1(1)=edge1(1)/2;
                end
                sqrt_=sqrt(edge3(2)^2+edge3(3)^2);
                if sqrt_>2*x_eps
                    x0=x0+x_eps*edge3/sqrt_;
                    edge3=edge3-2*x_eps*edge3/sqrt_;
                else
                    x0=x0+edge3/4;
                    edge3=edge3/2;
                end
                x0_=x0;
                while x0(2:3)*[e21;e31]<x_max || ...
                        (x0(2:3)+edge3(2:3))*[e21;e31]<x_max
                    h_=fill_box(x0,edge1,edge2,edge3,x_limit,...
                        pmt_display(m).color,pmt_display(m).alpha);
                    if h_~=0
                        h(m)=h_;
                    end
                    x0(2:3)=x0(2:3)+[ds21,ds31]*[e21;e31]*[e21,e31];
                end
                x0(2:3)=x0_(2:3)-[ds21,ds31]*[e21;e31]*[e21,e31];
                while x0(2:3)*[e21;e31]>-x_max || ...
                        (x0(2:3)+edge3(2:3))*[e21;e31]>-x_max
                    h_=fill_box(x0,edge1,edge2,edge3,x_limit,...
                        pmt_display(m).color,pmt_display(m).alpha);
                    if h_~=0
                        h(m)=h_;
                    end
                    x0(2:3)=x0(2:3)-[ds21,ds31]*[e21;e31]*[e21,e31];
                end
            end
        else
            % inhomogeneous stripe
            % Loop through structural blocks.
            c2=param(stripe.block{end}.c2,param_index)-1;
            for l3=1:length(stripe.block)
                block=stripe.block{l3};
                c2_=c2;
                c2=param(block.c2,param_index);
                m=block.pmt_index;
                if ~isempty(pmt_display(m).color)
                    x0=[x1_bot,c1*[ds21,ds31]+c2*[ds22,ds32]+[dx2,dx3]];
                    edge1=[x1_top-x1_bot,0,0];
                    edge2=[0,(c2_-c2)*[ds22,ds32]];
                    edge3=[0,(c1_-c1)*[ds21,ds31]*[e21;e31]*[e21,e31]];
                    if edge1(1)>2*x_eps
                        x0(1)=x0(1)+x_eps;
                        edge1(1)=edge1(1)-2*x_eps;
                    else
                        x0(1)=x0(1)+edge1(1)/4;
                        edge1(1)=edge1(1)/2;
                    end
                    sqrt_=sqrt(edge2(2)^2+edge2(3)^2);
                    if sqrt_>2*x_eps
                        x0=x0+x_eps*edge2/sqrt_;
                        edge2=edge2-2*x_eps*edge2/sqrt_;
                    else
                        x0=x0+edge2/4;
                        edge2=edge2/2;
                    end
                    sqrt_=sqrt(edge3(2)^2+edge3(3)^2);
                    if sqrt_>2*x_eps
                        x0=x0+x_eps*edge3/sqrt_;
                        edge3=edge3-2*x_eps*edge3/sqrt_;
                    else
                        x0=x0+edge3/4;
                        edge3=edge3/2;
                    end
                    x0_=x0;
                    while x0(2:3)*[e21;e31]<x_max || ...
                            (x0(2:3)+edge3(2:3))*[e21;e31]<x_max
                        x0__=x0;
                        while x0(2:3)*[e22;e32]<x_max || ...
                                (x0(2:3)+edge2(2:3))*[e22;e32]<x_max
                            h_=fill_box(x0,edge1,edge2,edge3,x_limit,...
                                pmt_display(m).color,pmt_display(m).alpha);
                            if h_~=0
                                h(m)=h_;
                            end
                            x0=x0+[0,ds22,ds32];
                        end
                        x0(2:3)=x0__(2:3)-[ds22,ds32];
                        while x0(2:3)*[e22;e32]>-x_max || ...
                                (x0(2:3)+edge2(2:3))*[e22;e32]>-x_max
                            h_=fill_box(x0,edge1,edge2,edge3,x_limit,...
                                pmt_display(m).color,pmt_display(m).alpha);
                            if h_~=0
                                h(m)=h_;
                            end
                            x0=x0-[0,ds22,ds32];
                        end
                        x0=x0__+[0,ds21,ds31];
                    end
                    x0(2:3)=x0_(2:3)-[ds21,ds31];
                    while x0(2:3)*[e21;e31]>-x_max || ...
                            (x0(2:3)+edge3(2:3))*[e21;e31]>-x_max
                        x0__=x0;
                        while x0(2:3)*[e22;e32]<x_max || ...
                                (x0(2:3)+edge2(2:3))*[e22;e32]<x_max
                            h_=fill_box(x0,edge1,edge2,edge3,x_limit,...
                                pmt_display(m).color,pmt_display(m).alpha);
                            if h_~=0
                                h(m)=h_;
                            end
                            x0=x0+[0,ds22,ds32];
                        end
                        x0(2:3)=x0__(2:3)-[ds22,ds32];
                        while x0(2:3)*[e22;e32]>-x_max || ...
                                (x0(2:3)+edge2(2:3))*[e22;e32]>-x_max
                            h_=fill_box(x0,edge1,edge2,edge3,x_limit,...
                                pmt_display(m).color,pmt_display(m).alpha);
                            if h_~=0
                                h(m)=h_;
                            end
                            x0=x0-[0,ds22,ds32];
                        end
                        x0=x0__-[0,ds21,ds31];
                    end
                end
            end
        end
    end
elseif stratum.type==3
    % coordinate break
    dx2=dx2+param(stratum.dx2,param_index);
    dx3=dx3+param(stratum.dx3,param_index);
else % stratum.type==4
    for count=1:stratum.rep_count
        for l1=1:length(stratum.stratum)
            [x1_top,dx2,dx3,h]=plot_stratum(...
                x1_top,dx2,dx3,h,stratum.stratum{l1},...
                param_index,pmt_display,x_limit,x_eps,x_max,...
                d21,d31,d22,d32);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function p=param(p,i)
% Index into multi-dimensional parameter. (Repmat extension of singleton
% dimensions is implicit.)
s=size(p);
n=length(s);
i(n+1:end)=[];
i(end+1:n)=1;
i(s==1)=1;
i=num2cell(i);
p=p(i{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function h=fill_box(x0,edge1,edge2,edge3,x_limit,c,a)
% Display a 3-D rectangle with corner vertex x0(1,1:3) and edge vectors
% edge1(1,1:3), edge2(1,1:3), and edge3(1,1:3). The rectangle is clipped to
% limits x_limit(1,1:3) (low limit) and x_limit(2,1:3) (high limit). The
% rgb display color is c(1,1:3), and opacity (alpha) is a. Return a handle
% to one of the faces, or 0 if all faces are entirely outside the limits.
h=0;
h_=fill_face([...
    x0;...
    x0+edge1;...
    x0+edge1+edge2;...
    x0+edge2...
    ],x_limit,c,a);
if h_~=0
    h=h_;
end
h_=fill_face([...
    x0+edge3;...
    x0+edge3+edge1;...
    x0+edge3+edge1+edge2;...
    x0+edge3+edge2...
    ],x_limit,c,a);
if h_~=0
    h=h_;
end
h_=fill_face([...
    x0;...
    x0+edge2;...
    x0+edge2+edge3;...
    x0+edge3...
    ],x_limit,c,a);
if h_~=0
    h=h_;
end
h_=fill_face([...
    x0+edge1;...
    x0+edge1+edge2;...
    x0+edge1+edge2+edge3;...
    x0+edge1+edge3...
    ],x_limit,c,a);
if h_~=0
    h=h_;
end
h_=fill_face([...
    x0;...
    x0+edge3;...
    x0+edge3+edge1;...
    x0+edge1...
    ],x_limit,c,a);
if h_~=0
    h=h_;
end
h_=fill_face([...
    x0+edge2;...
    x0+edge2+edge3;...
    x0+edge2+edge3+edge1;...
    x0+edge2+edge1...
    ],x_limit,c,a);
if h_~=0
    h=h_;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function h=fill_face(x,x_limit,c,a)
% Display a rectangle face with corner vertices x(:,1:3).
h=0;
x=clip_x(x,x_limit(1,1));
if isempty(x)
    return
end
x=-clip_x(-x,-x_limit(2,1));
if isempty(x)
    return
end
x=x(:,[2,3,1]);
x=clip_x(x,x_limit(1,2));
if isempty(x)
    return
end
x=-clip_x(-x,-x_limit(2,2));
if isempty(x)
    return
end
x=x(:,[2,3,1]);
x=clip_x(x,x_limit(1,3));
if isempty(x)
    return
end
x=-clip_x(-x,-x_limit(2,3));
if isempty(x)
    return
end
x=x(:,[2,3,1]);
h=patch(x(:,2),x(:,3),x(:,1),c,'FaceAlpha',a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function x=clip_x(x,x1_limit)
% Clip rectance face x(:,1:3) to limit x(:,1)>=x1_limit.
if all(x(:,1)>=x1_limit)
    return
end
if all(x(:,1)<=x1_limit)
    x=[];
    return
end
nx=size(x,1);
x_=[];
xa=x(end,:);
xb=x(1,:);
for j=1:nx
    xc=x(mod(j,nx)+1,:);
    if xb(1)>=x1_limit
        if xa(1)<x1_limit
            x_(end+1,:)=xa+(xb-xa)*((x1_limit-xa(1))/(xb(1)-xa(1)));
        end
        x_(end+1,:)=xb;
        if xc(1)<x1_limit
            x_(end+1,:)=xc+(xb-xc)*((x1_limit-xc(1))/(xb(1)-xc(1)));
        end
    end
    xa=xb;
    xb=xc;
end
x=x_;
