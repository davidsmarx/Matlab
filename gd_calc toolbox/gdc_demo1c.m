% gdc_demo1c: Uniperiodic, sinusoidal grating
%
% Test case from J. Opt. Soc. Am. A/Vol. 10, No. 12/Dec. 1993, pp.
% 2581-2591, Table 4 with h/d=100. (This script is adapted from
% gdc_demo1a.) Toggle the "if 1" branch for the TM case. The TE data in
% Table 4 corresponds to eff1, and TM corresponds to eff2.
%
%
% This demo program was tested with MATLAB R14 (SP2) on a Dell Precision
% 530 (vintage Feb 2002) with dual 2.2 GHz Xeon P4 processors and 2 Gb RAM.
%
%
% Output with "if 1" branch:
%
% >> gdc_demo1c
% Elapsed time is 28.734000 seconds.
%  
% Diffraction efficiencies (m1, eff1, eff2, eff3, eff4)
% R:
% -2    0.055787    0.008663    0.032225    0.032225
% -1    0.023269   0.0025455    0.012907    0.012907
%  0    0.022471   0.0017678    0.012119    0.012119
%
% Output with "if 0" branch:
%
% >> gdc_demo1c
% Elapsed time is 16.500000 seconds.
%  
% Diffraction efficiencies (m1, eff1, eff2, eff3, eff4)
% R:
% -2     0.18472   0.0057095    0.095215    0.095215
% -1    0.027206    0.016335     0.02177     0.02177
%  0      0.0158    0.074683    0.045241    0.045241
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

demo=true; % true: run gdc_demo_engine; false: run gdc_engine

% If ctr_sect==true, define stratum widths by center-sectioning; otherwise
% use (slightly) more accurate averaging.
ctr_sect=true;

% Define parameters for grating structure and incident field:
grating_pmt=(0.3+7.0*i)^2; % grating permittivity
d=1.0; % grating period
wavelength=d/1.7;
h=100*d; % grating height
if 1
    % TE case in JOSA publication
    L1=50; % number of grating strata (for "staircase" approximation)
    m_max=25; % max diffraction order index
else
    % TM case in JOSA publication
    L1=10; % number of grating strata (for "staircase" approximation)
    m_max=52; % max diffraction order index
end
% theta1 = angle between the grating normal (x1 axis) and the incident wave
% vector
theta1=30*pi/180;
% phi1 = angle between the grating-tangential direction normal to the
% grating lines (x2 axis) and the incident wave vector's grating-tangential
% projection (x2, x3 plane)
phi1=0*pi/180;

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
% projection is implicitly f1=-cos(theta1)/wavelength.
inc_field.f2=sin(theta1)*cos(phi1)/wavelength;
inc_field.f3=sin(theta1)*sin(phi1)/wavelength;

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

% Compute the diffraction efficiencies. (Only the reflected waves are
% retained.)
R=gdc_eff(scat_field,inc_field);
% Discard diffracted waves that decay exponentially with distance from the
% grating. (These include evanescent waves.)
% (Note: "[scat_field.f1r]" is MATLAB syntax for
% "scat_field(1).f1r, scat_field(2).f1r, ...".)
R=R(imag([scat_field.f1r])==0);
% Extract the diffraction order indices for the reflected waves. (Only the
% m1 indices matter; the m2's are all zero.)
R_m1=[R.m1].';
% Extract the reflection efficiencies (R1...R4) corresponding to the four
% incident polarization states defined in gdc_eff.m.
R1=[R.eff1].';
R2=[R.eff2].';
R3=[R.eff3].';
R4=[R.eff4].';
% Tabulate the diffraction order indices and diffraction efficiencies.
disp(' ');
disp('Diffraction efficiencies (m1, eff1, eff2, eff3, eff4)');
disp('R:');
disp(num2str([R_m1 R1 R2 R3 R4]));
