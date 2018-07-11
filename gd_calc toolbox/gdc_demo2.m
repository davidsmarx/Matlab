% gdc_demo2: Biperiodic grating comprising rectangular pyramids
%
% Test case from G. H. Derrick et al, Applied Physics 18, 39-52 (1979).
% See also J. J. Greffet et al, Optics Letters 17(24), 1740-1743 (1992);
% R. Brauer and O. Bryngdahl, Optics Communications 100, 1-5 (1993).
%
% This demo covers the following topics:
%   - GDC interface and setup for a simple biperiodic grating
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
% >> gdc_demo2
% Elapsed time is 31.827276 seconds.
%  
% Diffraction efficiencies (m1, m2, eff2)
% R:
% -1           0   0.0024859
%  0           0    0.019436
% T:
% -1          -1  0.00085096
%  0          -1   0.0067344
% -1           0   0.0029004
%  0           0     0.96485
%  1           0   0.0027419
% Energy loss:
% -1.3856e-013
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

demo=true; % true: run gdc_demo_engine; false: run gdc_engine

% Define parameters for grating structure and incident field:
d1=1.50; % grating period in x2 direction
d2=1.00; % grating period in x3 direction
h=0.25; % grating height
wavelength=1.533;
phi=30*pi/180; % incidence polar angle from x1 axis
psi=45*pi/180; % incidence azimuthal angle around x1 axis (zero on x2 axis)
grating_pmt=1.5^2; % grating permittivity
L1=16; % number of grating strata
m_max=5; % maximum diffraction order index

if ~demo % The following code block is replicated in gdc_demo_engine.
    % Construct grating.
    clear grating
    grating.pmt={1,grating_pmt}; % grating material permittivities
    grating.pmt_sub_index=2; % substrate permittivity index
    grating.pmt_sup_index=1; % superstrate permittivity index
    grating.d21=d1; % first grating period: x2 projection
    grating.d31=0; % first grating period: x3 projection
    grating.d22=0; % second grating period: x2 projection
    grating.d32=d2; % second grating period: x3 projection
    grating.stratum={};
    clear stratum
    stratum.type=2; % biperiodic stratum
    stratum.thick=h/L1; % stratum thickness
    % The following h11 ... h22 spec indicates that the stratum's period
    % vectors match the grating periods (GD-Calc.pdf, equation 3.18).
    stratum.h11=1;
    stratum.h12=0;
    stratum.h21=0;
    stratum.h22=1;
    clear stripe
    stripe.type=0; % first stripe is homogeneous
    stripe.pmt_index=1; % first stripe's permittivity index
    stratum.stripe{1}=stripe;
    clear stripe
    stripe.type=1; % second stripe is inhomogeneous 
    stripe.block{1}.pmt_index=1; % first block's permittivity index
    stripe.block{2}.pmt_index=2; % second block's permittivity index
    stratum.stripe{2}=stripe;
    clear stripe
    for l1=1:L1
        % Set first and second stripes' boundaries (c1) on positive side:
        stratum.stripe{1}.c1=-((L1-l1+0.5)/L1)/2;
        stratum.stripe{2}.c1=-stratum.stripe{1}.c1;
        % Set first and second block' boundaries (c2) on positive side:
        stratum.stripe{2}.block{1}.c2=-((L1-l1+0.5)/L1)/2;
        stratum.stripe{2}.block{2}.c2=-stratum.stripe{2}.block{1}.c2;
        grating.stratum{end+1}=stratum;
    end
    clear stratum
end

% Define the indicent field.
clear inc_field
inc_field.wavelength=wavelength;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector. The grating-normal (x1)
% projection is implicitly f1=-cos(phi)/wavelength.
inc_field.wavelength=wavelength;
inc_field.f2=sin(phi)*cos(psi)/wavelength;
inc_field.f3=sin(phi)*sin(psi)/wavelength;

% Specify which diffracted orders are to be retained in the calculations.
order=[];
for m2=-m_max:m_max
    order(end+1).m2=m2;
    order(end).m1=-m_max:m_max;
end

% Run the diffraction calculations.
tic
if demo
    [grating,param_size,scat_field,inc_field]=...
        gdc_demo_engine(2,inc_field,order,...
        grating_pmt,d1,d2,h,L1);
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
% Tabulate the diffraction order indices and diffraction efficiencies for
% an incident field polarized parallel to the incident plane. Also tabulate
% the fractional energy loss in the grating.
disp(' ');
disp('Diffraction efficiencies (m1, m2, eff2)');
disp('R:');
disp(num2str([[R.m1].' [R.m2].' [R.eff2].']));
disp('T:');
disp(num2str([[T.m1].' [T.m2].' [T.eff2].']));
disp('Energy loss:');
disp(num2str(1-sum([[[R.eff2].']; [[T.eff2].']])));

% Plot the grating.
clear pmt_display
pmt_display(1).name='';
pmt_display(1).color=[];
pmt_display(1).alpha=1;
pmt_display(2).name='';
pmt_display(2).color=[.75,.75,.75];
pmt_display(2).alpha=1;
x_limit=[-0.5*h,-1.5*d1,-1.5*d2;1.5*h,1.5*d1,1.5*d2];
h_plot=gdc_plot(grating,1,pmt_display,x_limit);
if ~isempty(h_plot)
    view(-25,15)
end
