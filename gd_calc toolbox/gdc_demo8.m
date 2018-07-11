% gdc_demo8: Biperiodic grating - square metal grid
%
% Test case from J. Opt. A: Pure Appl. Opt. 6 (2004), pp. 43-50, Example 3.
% (The published diffractions efficiencies correspond to eff3 in the output
% listing below.)
%
% gdc_demo8 covers the following topics:
%	- biperiodic metallic grating
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
% >> gdc_demo8
% Elapsed time is 759.334817 seconds.
%  
% Diffraction efficiencies (m1, m2, eff1, eff2, eff3, eff4)
% R:
%  0          -1    0.027187    0.044984    0.036085    0.036085
% -1           0       0.045    0.027216    0.036108    0.036108
%  0           0     0.13873     0.13872     0.13872     0.13872
%  1           0       0.045    0.027216    0.036108    0.036108
%  0           1    0.027187    0.044984    0.036085    0.036085
% T:
%  0          -1   0.0071299    0.075834    0.041482    0.041482
% -1           0    0.075847   0.0071583    0.041503    0.041503
%  0           0     0.44691     0.44673     0.44682     0.44682
%  1           0    0.075847   0.0071583    0.041503    0.041503
%  0           1   0.0071299    0.075834    0.041482    0.041482
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

demo=true; % true: run gdc_demo_engine; false: run gdc_engine

% Define parameters for grating structure and incident field:
stratum_pmt=(1+5*i)^2; % metal grid permittivity
wavelength=1;
h=0.2*wavelength; % grating height
d=1.2*wavelength; % grating period
c=5/6; % hole width/period
m_max=13; % maximum diffraction order index

if ~demo % The following code block is replicated in gdc_demo_engine.
    % Construct grating.
    clear grating
    grating.pmt={1,stratum_pmt}; % grating material permittivities
    grating.pmt_sub_index=1; % substrate permittivity index
    grating.pmt_sup_index=1; % superstrate permittivity index
    grating.d21=d; % first grating period: x2 projection
    grating.d31=0; % first grating period: x3 projection
    grating.d22=0; % second grating period: x2 projection
    grating.d32=d; % second grating period: x3 projection
    % Construct the stratum.
    clear stratum
    stratum.type=2; % biperiodic stratum
    stratum.thick=h; % stratum thickness
    % The following h11 ... h22 spec indicates that the stratum's period
    % vectors match the grating periods (GD-Calc.pdf, equation 3.18).
    stratum.h11=1;
    stratum.h12=0;
    stratum.h21=0;
    stratum.h22=1;
    stratum.stripe={};
    clear stripe block
    stripe.type=0; % first stripe is homogeneous
    stripe.c1=-c/2; % first stripe's boundary on positive side
    stripe.pmt_index=2; % first stripe's permittivity index
    stratum.stripe{end+1}=stripe;
    clear stripe
    stripe.type=1; % second stripe is homogeneous
    stripe.c1=c/2; % second stripe's boundary on positive side
    block.c2=-c/2; % first block's boundary on positive side
    block.pmt_index=2; % first block's permittivity index
    stripe.block{1}=block;
    block.c2=c/2; % second block's boundary on positive side
    block.pmt_index=1; % second block's permittivity index
    stripe.block{2}=block;
    clear block
    stratum.stripe{end+1}=stripe;
    clear stripe
    grating.stratum={stratum};
    clear stratum
end

% Define the indicent field (normal incidence).
clear inc_field
inc_field.wavelength=wavelength;
inc_field.f2=0;
inc_field.f3=0;

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
        gdc_demo_engine(8,inc_field,order,...
        stratum_pmt,d,h,c);
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
% the four incident polarization states defined in gdc_eff.m.
disp(' ');
disp('Diffraction efficiencies (m1, m2, eff1, eff2, eff3, eff4)');
disp('R:');
disp(num2str([[R.m1].' [R.m2].' ...
    [R.eff1].' [R.eff2].' [R.eff3].' [R.eff4].']));
disp('T:');
disp(num2str([[T.m1].' [T.m2].' ...
    [T.eff1].' [T.eff2].' [T.eff3].' [T.eff4].']));

% Plot the grating.
clear pmt_display
pmt_display(1).name='';
pmt_display(1).color=[];
pmt_display(1).alpha=1;
pmt_display(2).name='';
pmt_display(2).color=[.75,.75,.75];
pmt_display(2).alpha=1;
x_limit=[-0.5*h,-1.5*d,-1.5*d;1.5*h,1.5*d,1.5*d];
gdc_plot(grating,1,pmt_display,x_limit);
