% gdc_demo7: Biperiodic grating - skewed metal grid, unit cell B
%
% Test case from J. Opt. Soc. Am. A/Vol. 14, No. 10/Oct. 1997, pp.
% 2758-2767, Example 3, but with a different unit cell.
%
% gdc_demo6 and gdc_demo7 cover the following topics:
%   - biperiodic grating modeling with different unit cells and stripe
%     orientations
%   - biperiodic grating with non-orthogonal period vectors
%	- metallic grating
%   - diffraction order selection
%   - Compare computation results with published numeric data.
%     (Note: The results do not agree with the published data, but the
%     published reflectance efficiencies in the JOSA paper appear
%     unexpectedly high in comparison to a bare metal substrate
%     reflectivity.)
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf (Part 1)
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
%
% This demo program was tested with MATLAB R14 (SP2) on a Dell Precision
% 530 (vintage Feb 2002) with dual 2.2 GHz Xeon P4 processors and 2 Gb RAM.
%
% Output:
%
% >> gdc_demo7
% Elapsed time is 738.391000 seconds.
%  
% Diffraction efficiencies (m1, m2, eff2)
% R:
% -2          -2  0.00037692
% -1          -2  0.00059872
% -2          -1  0.00059872
% -1          -1   0.0050594
%  0          -1    0.009955
% -1           0    0.009955
%  0           0     0.11032
%  1           0   0.0063215
%  0           1   0.0063215
% T:
% -3          -3  1.349e-005
% -2          -3  0.00014536
% -1          -3  0.00027968
% -3          -2  0.00014536
% -2          -2  0.00052813
% -1          -2   0.0014759
%  0          -2  0.00098662
% -3          -1  0.00027968
% -2          -1   0.0014759
% -1          -1   0.0039865
%  0          -1   0.0032974
%  1          -1   0.0006867
% -2           0  0.00098662
% -1           0   0.0032974
%  0           0   0.0031037
%  1           0  0.00070066
%  2           0 9.0797e-005
% -1           1   0.0006867
%  0           1  0.00070066
%  1           1  0.00014042
%  2           1 2.0585e-005
%  0           2 9.0797e-005
%  1           2 2.0585e-005
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

demo=true; % true: run gdc_demo_engine; false: run gdc_engine

alt_stripe=1; % selector for alternate stripe orientation (1 or 2)

% Define parameters for grating structure and incident field:
substrate_pmt=2.25; % substrate permittivity
stratum_pmt=1+5*i; % metal grid permittivity
wavelength=1;
h=wavelength; % grating height
d=2*wavelength; % grating period
theta=30*pi/180; % incidence polar angle from x1 axis
phi=30*pi/180; % incidence azimuthal angle around x1 axis (zero on x2 axis)
N=1000; % number of "teeth" on each edge of the hole = N-1
m_max=12; % maximum diffraction order index

if ~demo % The following code block is replicated in gdc_demo_engine.
    % Construct grating.
    clear grating stratum stripe block
    grating.pmt=...
        {1,substrate_pmt,stratum_pmt}; % grating material permittivities
    grating.pmt_sub_index=2; % substrate permittivity index
    grating.pmt_sup_index=1; % superstrate permittivity index
    % Define the x2 and x3 projections of the first grating period
    % (d21,d31) and second grating period (d22,d32), corresponding to unit
    % cell B. The second period is parallel to the grid holes' long
    % diagonals, and the first period is parallel to the short diagonals.
    % (Note: The second period defines the stripe orientation - see Figure
    % 3 in GD-Calc.pdf.)
    zeta=30*pi/180;
    grating.d21=-d*sin(zeta);
    grating.d31=d*cos(zeta);
    grating.d22=d*(1+sin(zeta));
    grating.d32=d*cos(zeta);
    % Construct the stratum. (Refer to Figure 3 in GD-Calc.pdf and view the
    % grating plot with a small N value to follow the construction logic.)
    % The x2, x3 coordinate origin is at the center of a grid hole.
    stratum.type=2;
    stratum.thick=h;
    if alt_stripe==1
        % Stripes are parallel to [grating.d22,grating.d32].
        stratum.h11=1;
        stratum.h12=0;
        stratum.h21=0;
        stratum.h22=1;
    else
        % Stripes are parallel to [grating.d21,grating.d31].
        stratum.h11=0;
        stratum.h12=1;
        stratum.h21=1;
        stratum.h22=0;
    end
    stratum.stripe=cell(1,4*N-2);
    stripe.type=1;
    dc1=0.25/(N-0.5);
    dc2=dc1;
    stripe.c1=-0.25+dc1;
    block.c2=-dc2/2;
    block.pmt_index=3;
    stripe.block{1}=block;
    block.c2=-block.c2;
    block.pmt_index=1;
    stripe.block{2}=block;
    for l2=1:N-1
        stratum.stripe{l2}=stripe;
        stripe.c1=stripe.c1+dc1;
        stripe.block{1}.c2=stripe.block{1}.c2-dc2;
        stripe.block{2}.c2=-stripe.block{1}.c2;
    end
    stripe.block{1}.c2=stripe.block{1}.c2+dc2/4;
    stripe.block{2}.c2=-stripe.block{1}.c2;
    stratum.stripe{N}=stripe;
    for l2=N+1:2*N-1
        stripe=stratum.stripe{2*N-l2};
        stripe.c1=-stripe.c1+dc1;
        stratum.stripe{l2}=stripe;
    end
    for l2=2*N:4*N-2
        stripe=stratum.stripe{l2-2*N+1};
        stripe.c1=stripe.c1+0.5;
        stripe.block{1}.c2=stripe.block{1}.c2+0.5;
        stripe.block{2}.c2=stripe.block{2}.c2+0.5;
        stratum.stripe{l2}=stripe;
    end
    grating.stratum={stratum};
    clear stratum stripe block
end

% Define the indicent field.
clear inc_field
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector. The grating-normal (x1)
% projection is implicitly f1=-cos(theta)/wavelength.
inc_field.wavelength=wavelength;
inc_field.f2=sin(theta)*cos(phi)/wavelength;
inc_field.f3=sin(theta)*sin(phi)/wavelength;

% Specify which diffracted orders are to be retained in the calculations.
% First define order indices relative to the unit cell A.
m1=repmat(-m_max:m_max,[2*m_max+1,1]);
m2=m1.';
m1=m1(:).';
m2=m2(:).';
% Redefine order indices relative to unit cell B.
[m1,m2]=deal(m1-m2,m1+m2);
% Construct the order struct.
order=[];
for m2_=unique(m2)
    order(end+1).m2=m2_;
    order(end).m1=m1(m2==m2_);
end
clear m1 m2

% Run the diffraction calculations.
tic
if demo
    [grating,param_size,scat_field,inc_field]=...
        gdc_demo_engine(7,inc_field,order,...
        substrate_pmt,stratum_pmt,d,h,N,alt_stripe);
else
    [param_size,scat_field,inc_field]=gdc(grating,inc_field,order);
end
toc
if isempty(scat_field)
    disp('Interrupted by user.');
    return
end

% Compute the diffraction efficiencies.
[R,T]=gdc_eff(scat_field,inc_field);
% Discard diffracted waves that decay exponentially with distance from the
% grating. (These include evanescent waves and, if the substrate's
% permittivity is not real-valued, all transmitted waves.)
R=R(imag([scat_field.f1r])==0);
T=T(imag([scat_field.f1t])==0);
% Extract the diffraction order indices for the reflected and transmitted
% waves.
R_m1=[R.m1].';
R_m2=[R.m2].';
T_m1=[T.m1].';
T_m2=[T.m2].';
% Redefine order indices relative to unit cell B.
[R_m1,R_m2]=deal((R_m1+R_m2)/2,(R_m2-R_m1)/2);
[T_m1,T_m2]=deal((T_m1+T_m2)/2,(T_m2-T_m1)/2);
% Sort the orders by order index
[R_m1,k]=sort(R_m1);
R_m2=R_m2(k);
R=R(k);
[R_m2,k]=sort(R_m2);
R_m1=R_m1(k);
R=R(k);
[T_m1,k]=sort(T_m1);
T_m2=T_m2(k);
T=T(k);
[T_m2,k]=sort(T_m2);
T_m1=T_m1(k);
T=T(k);
% Tabulate the diffraction order indices and diffraction efficiencies for
% an incident field polarized parallel to the incident plane.
disp(' ');
disp('Diffraction efficiencies (m1, m2, eff2)');
disp('R:');
disp(num2str([R_m1 R_m2 [R.eff2].']));
disp('T:');
disp(num2str([T_m1 T_m2 [T.eff2].']));

% Plot the grating. (Caution: Cancel the plot if N is very large.)
clear pmt_display
pmt_display(1).name='';
pmt_display(1).color=[];
pmt_display(1).alpha=1;
pmt_display(2).name='';
pmt_display(2).color=[.5,.5,.5];
pmt_display(2).alpha=1;
pmt_display(3).name='';
pmt_display(3).color=[.75,.75,.75];
pmt_display(3).alpha=1;
x_limit=[-0.5*h,-1.5*d,-1.5*d;1.5*h,1.5*d,1.5*d];
h_plot=gdc_plot(grating,1,pmt_display,x_limit);
if ~isempty(h_plot)
    view(-45,45)
end
