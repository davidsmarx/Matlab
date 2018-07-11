% gdc_demo4: Biperiodic checkerboard grating
%
% Test case from J. Opt. Soc. Am. A/Vol. 14, No. 10/Oct. 1997, pp.
% 2758-2767, Example 1 with unit cell A.
%
% gc_demo4 covers the following topics:
%   - GDC interface and setup for a simple biperiodic grating
%   - choice of unit cell and stripe orientations
%   - diffraction order selection
%   - stripe-induced symmetry error
%   - Compare computation results with published numeric data.
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
% >> gdc_demo4
% Elapsed time is 127.750000 seconds.
%
% Diffraction efficiencies (with incident E parallel to e2 or e3)
% R:
% -1          -1   0.0012662   0.0005683
%  0          -1   0.0016871   0.0016883
%  1          -1  0.00056848   0.0012629
% -1           0   0.0016871   0.0016883
%  0           0    0.011109    0.011111
%  1           0   0.0016871   0.0016883
% -1           1  0.00056848   0.0012629
%  0           1   0.0016871   0.0016883
%  1           1   0.0012662   0.0005683
% T:
% -1          -2   0.0096013   0.0061868
%  0          -2   0.0067862   0.0067372
%  1          -2   0.0062164    0.009576
% -2          -1   0.0096013   0.0061868
% -1          -1    0.040219    0.054862
%  0          -1     0.13067      0.1306
%  1          -1    0.054803    0.040345
%  2          -1   0.0062164    0.009576
% -2           0   0.0067862   0.0067372
% -1           0     0.13067      0.1306
%  0           0     0.17532     0.17565
%  1           0     0.13067      0.1306
%  2           0   0.0067862   0.0067372
% -2           1   0.0062164    0.009576
% -1           1    0.054803    0.040345
%  0           1     0.13067      0.1306
%  1           1    0.040219    0.054862
%  2           1   0.0096013   0.0061868
% -1           2   0.0062164    0.009576
%  0           2   0.0067862   0.0067372
%  1           2   0.0096013   0.0061868
% Energy loss:
% -4.0612e-013 3.9635e-013
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

demo=true; % true: run gdc_demo_engine; false: run gdc_engine

alt_order=false; % Toggle for alternate order truncation.

% Define parameters for grating structure and incident field:
grating_pmt=2.25; % grating permittivity
wavelength=1;
h=wavelength; % grating height
w=1.25*wavelength; % width of squares (half-period along x2, x3 axes)
m_max=10; % maximum diffraction order index

if ~demo % The following code block is replicated in gdc_demo_engine.
    % Construct grating.
    clear grating
    grating.pmt={1,grating_pmt}; % grating material permittivities
    grating.pmt_sub_index=2; % substrate permittivity index
    grating.pmt_sup_index=1; % superstrate permittivity index
    % Define the x2 and x3 projections of the first grating period
    % (d21,d31) and second grating period (d22,d32). The x2 and x3 axes are
    % aligned to unit cell A, and the period vectors define the unit cell's
    % edges.
    grating.d21=2*w;
    grating.d31=0;
    grating.d22=0;
    grating.d32=2*w;
    % Construct the stratum. (View the grating plot plot to follow the
    % construction logic.)
    clear stratum stripe block
    stratum.type=2;
    stratum.thick=h;
    stratum.h11=1;
    stratum.h12=0;
    stratum.h21=0;
    stratum.h22=1;
    stripe.c1=0.25;
    stripe.type=1;
    block.c2=0.25;
    block.pmt_index=2;
    stripe.block{1}=block;
    block.c2=0.75;
    block.pmt_index=1;
    stripe.block{2}=block;
    stratum.stripe{1}=stripe;
    stripe.c1=0.75;
    block.c2=0.25;
    stripe.block{1}=block;
    block.c2=0.75;
    block.pmt_index=2;
    stripe.block{2}=block;
    stratum.stripe{2}=stripe;
    grating.stratum{1}=stratum;
    clear stratum stripe block
end

% Define the indicent field (normal incidence).
clear inc_field
inc_field.wavelength=wavelength;
inc_field.f2=0;
inc_field.f3=0;

% Specify which diffracted orders are to be retained in the calculations.
% First define order indices relative to unit cell B.
m1=repmat(-m_max:m_max,[2*m_max+1,1]);
m2=m1.';
m1=m1(:).';
m2=m2(:).';
if alt_order
    % Alternate order truncation: |m1|+|m2|<=m_max
    find_=find(abs(m1)+abs(m2)<=m_max);
    m1=m1(find_);
    m2=m2(find_);
end
% Redefine order indices relative to unit cell A.
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
        gdc_demo_engine(4,inc_field,order,...
        grating_pmt,w,h);
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
% Extract the reflection efficiencies (R_e2,R_e3) and transmission
% efficiencies (T_e2,T_e3), wherein R_e2 and T_e2 correspond to an incident
% field polarized in the e2 direction, and R_e3 and T_e3 correspond to an
% incident field polarized in the e3 direction. For the current diffraction
% geometry, e2 and e3 are equal to the polarization basis vectors p and s,
% respectively, so the e2 and e3 polarization efficiencies correspond
% respectively to eff2 and eff1.
R_e2=[R.eff2].';
R_e3=[R.eff1].';
T_e2=[T.eff2].';
T_e3=[T.eff1].';

% Tabulate the diffraction order indices and diffraction efficiencies for
% an incident field polarized parallel to the x2 or x3 axis (unit basis
% vector e2 or e3). Also tabulate the fractional energy loss in the
% grating. (Note: For the present diffraction geometry the incident field's
% polarization basis vectors s and p are aligned to e3 and e2,
% respectively - see GD-Calc.pdf, equations 4.18 and 4.19.)
disp(' ');
disp('Diffraction efficiencies (with incident E parallel to e2 or e3)');
disp('R:');
disp(num2str([R_m1 R_m2 R_e2 R_e3]));
disp('T:');
disp(num2str([T_m1 T_m2 T_e2 T_e3]));
disp('Energy loss:');
disp(num2str(1-sum([[R_e3 R_e2]; [T_e3 T_e2]])));

% Plot the grating.
clear pmt_display
pmt_display(1).name='';
pmt_display(1).color=[];
pmt_display(1).alpha=1;
pmt_display(2).name='';
pmt_display(2).color=[.75,.75,.75];
pmt_display(2).alpha=1;
x_limit=[-0.5*h,-2.5*w,-2.5*w;1.5*h,2.5*w,2.5*w];
h_plot=gdc_plot(grating,1,pmt_display,x_limit);
if ~isempty(h_plot)
    view(-30,45)
end
