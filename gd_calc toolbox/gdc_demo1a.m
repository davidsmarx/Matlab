% gdc_demo1a: Uniperiodic, sinusoidal grating
%
% Test case from J. Opt. Soc. Am. A/Vol. 10, No. 12/Dec. 1993, pp.
% 2581-2591, Table 2. (N in Table 2 is 2*m_max+1, and M is L1.)
% (This demo can also be modified to match the grating configuration for
% Tables 3 and 4 - see gdc_demo1b and gdc_demo1c.)
%
% This demo covers the following topics:
%   - Illustrate GDC interface and setup for a simple uniperiodic grating.
%   - Compute conical diffraction with a non-trivial polarization geometry.
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
% >> gdc_demo1a
% Elapsed time is 13.344000 seconds.
%  
% Diffraction efficiencies (m1, eff1, eff2, eff3, eff4)
% R:
% -3    0.011769    0.020237    0.013569    0.016679
% -2    0.039664    0.033315    0.029423    0.036985
% -1    0.040045   0.0054773    0.020138    0.022118
%  0     0.10368    0.010056    0.061266    0.058194
% T:
% -5  0.00018774  0.00044337  0.00030678  0.00023839
% -4 5.1896e-005  0.00054232  0.00016879  0.00029671
% -3   0.0073647    0.016271    0.011378    0.011246
% -2    0.049306    0.097487    0.070731    0.070444
% -1    0.099548     0.17061     0.12809     0.12797
%  0    0.070222    0.048023    0.065278    0.057236
%  1     0.51523     0.55687     0.54436     0.54243
%  2     0.06294    0.040671    0.055285    0.056159
% Energy loss:
% 1.4433e-015 -5.3646e-013 -2.7978e-013 -2.4736e-013
%  
% Diffraction efficiencies (with H3=0, E3=0)
% R:
% -3    0.011242    0.020764
% -2    0.037459     0.03552
% -1    0.038523   0.0069991
%  0     0.10292    0.010815
% T:
% -5  0.00019076  0.00044035
% -4 2.5097e-005  0.00056912
% -3   0.0074306    0.016205
% -2    0.049579    0.097213
% -1    0.099059      0.1711
%  0    0.071536    0.046708
%  1     0.51857     0.55353
%  2     0.06347     0.04014
% Energy loss:
% -1.3989e-014 -5.2136e-013
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

demo=false; % true: run gdc_demo_engine; false: run gdc_engine

% If ctr_sect==true, define stratum widths by center-sectioning; otherwise
% use (slightly) more accurate averaging.
ctr_sect=true;

% Define parameters for grating structure and incident field:
grating_pmt=4.; % grating permittivity
d=1.0; % grating period
wavelength=d/2.0;
h=0.6*wavelength; % grating height
L1=50; % number of grating strata (for "staircase" approximation)
m_max=20; % max diffraction order index
% theta3 = angle between the grating lines (x3 axis) and the incident wave
% vector
theta3=75*pi/180;
% phi3 = angle between the grating normal (x1 axis) and the incident wave
% vector's projection normal to the grating lines (i.e., x1, x2 plane
% projection)
phi3=60*pi/180;

if ~demo % The following code block is replicated in gdc_demo_engine.
    if ctr_sect
        % Define stratum half-widths (in units of d) by center-sectioning.
        c1=acos(-1+2*((1:L1)-0.5)/L1)/(2*pi);
    else
        % Define stratum half-widths using slightly more accurate averaging
        % method.
        c1=2*(1:L1-1)/L1-1;
        c1=diff([-0.25, (c1.*acos(c1)-sqrt(1-c1.^2))/(4*pi), 0])*L1;
    end

    % Construct grating. (Note: Two grating periods must be specified
    % because the grating can generally be biperiodic, although for this
    % example the second period (d22,d32) is irrelevant.)
    clear grating
    grating.pmt={1,grating_pmt}; % grating material permittivities
    grating.pmt_sub_index=2; % substrate permittivity index
    grating.pmt_sup_index=1; % superstrate permittivity index
    grating.d21=d; % first grating period: x2 projection
    grating.d31=0; % first grating period: x3 projection
    grating.d22=0; % second grating period: x2 projection
    grating.d32=d; % second grating period: x3 projection
    grating.stratum={};
    clear stratum stripe
    stratum.type=1; % uniperiodic stratum
    stratum.thick=h/L1; % stratum thickness
    % The following h11, h12 spec indicates that the stratum's period
    % vector matches the first grating period (GD-Calc.pdf, equations 3.22
    % and 3.23).
    stratum.h11=1;
    stratum.h12=0;
    for l1=1:L1
        stripe.pmt_index=1; % first stripe's permittivity index
        stripe.c1=-c1(l1); % first stripe's boundary on positive side
        stratum.stripe{1}=stripe;
        stripe.pmt_index=2; % second stripe's permittivity index
        stripe.c1=c1(l1); % second stripe's boundary on positive side
        stratum.stripe{2}=stripe;
        grating.stratum{end+1}=stratum;
    end
    clear c1 stratum stripe
end

% Define the indicent field.
clear inc_field
inc_field.wavelength=wavelength;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector. The grating-normal (x1)
% projection is implicitly f1=-sin(theta3)*cos(phi3)/wavelength.
inc_field.f2=sin(theta3)*sin(phi3)/wavelength;
inc_field.f3=cos(theta3)/wavelength;

% Specify which diffracted orders are to be retained in the calculations.
% (m2 is zero because this is a uniperiodic grating - all diffraction
% orders for m2~=0 have zero amplitude.)
clear order
order(1).m2=0;
order(1).m1=-m_max:m_max;

% Run the diffraction calculations.
tic
if demo
    [grating,param_size,scat_field,inc_field]=...
        gdc_demo_engine(1,inc_field,order,...
        grating_pmt,d,h,L1,ctr_sect);
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
% (Note: "[scat_field.f1r]" is MATLAB syntax for
% "scat_field(1).f1r, scat_field(2).f1r, ...".)
R=R(imag([scat_field.f1r])==0);
T=T(imag([scat_field.f1t])==0);
% Extract the diffraction order indices for the reflected and transmitted
% waves. (Only the m1 indices matter; the m2's are all zero.)
R_m1=[R.m1].';
T_m1=[T.m1].';
% Extract the reflection efficiencies (R1...R4) and transmission
% efficiencies (T1...T4) corresponding to the four incident polarization
% states defined in gdc_eff.m.
R1=[R.eff1].';
R2=[R.eff2].';
R3=[R.eff3].';
R4=[R.eff4].';
T1=[T.eff1].';
T2=[T.eff2].';
T3=[T.eff3].';
T4=[T.eff4].';
% Tabulate the diffraction order indices and diffraction efficiencies. Also
% tabulate the fractional energy loss in the grating.
disp(' ');
disp('Diffraction efficiencies (m1, eff1, eff2, eff3, eff4)');
disp('R:');
disp(num2str([R_m1 R1 R2 R3 R4]));
disp('T:');
disp(num2str([T_m1 T1 T2 T3 T4]));
disp('Energy loss:');
disp(num2str(1-sum([[R1 R2 R3 R4]; [T1 T2 T3 T4]])));

% Compute the diffraction efficiencies for the incident E field's
% polarization normal to the grating lines (x3 axis). Assume an incident E
% field complex amplitude of the form A*s+B*p, as described in gdc_eff.m,
% with A and B defined so that the field is orthogonal to the x3 axis.
A=-inc_field.p3;
B=inc_field.s3;
% Normalize abs(A)^2+abs(B)^2 to 1.
C=1/sqrt(abs(A)^2+abs(B)^2);
A=A*C;
B=B*C;
% Calculate the reflection and transmission efficiencies, denoted R5 and
% T5, for the above-defined incident field.
R5=abs(A)^2*R1+abs(B)^2*R2...
    +real(conj(A)*B)*(2*R3-R1-R2)-imag(conj(A)*B)*(2*R4-R1-R2);
T5=abs(A)^2*T1+abs(B)^2*T2...
    +real(conj(A)*B)*(2*T3-T1-T2)-imag(conj(A)*B)*(2*T4-T1-T2);
% Do the same for the orthogonal incident polarization (H field orthogonal
% to the x3 axis); denote the corresponding efficiencies as R6, T6:
[A,B]=deal(B,-A);
R6=abs(A)^2*R1+abs(B)^2*R2...
    +real(conj(A)*B)*(2*R3-R1-R2)-imag(conj(A)*B)*(2*R4-R1-R2);
T6=abs(A)^2*T1+abs(B)^2*T2...
    +real(conj(A)*B)*(2*T3-T1-T2)-imag(conj(A)*B)*(2*T4-T1-T2);
% Tabulate the results.
disp(' ');
disp('Diffraction efficiencies (with H3=0, E3=0)');
disp('R:');
disp(num2str([R_m1 R6 R5]));
disp('T:');
disp(num2str([T_m1 T6 T5]));
disp('Energy loss:');
disp(num2str(1-sum([[R6 R5]; [T6 T5]])));

% Plot the grating.
clear pmt_display
pmt_display(1).name='';
pmt_display(1).color=[];
pmt_display(1).alpha=1;
pmt_display(2).name='';
pmt_display(2).color=[.75,.75,.75];
pmt_display(2).alpha=1;
x_limit=[-0.5*h,-d,-d;1.5*h,d,d];
h_plot=gdc_plot(grating,1,pmt_display,x_limit);
if ~isempty(h_plot)
    view(127.5,20) % Match view in GD-Calc_Demo.pdf.
end
