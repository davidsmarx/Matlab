% gdc_demo10: slanted lamellar grating
%
% This demo models a slanted lamellar grating configured to operate as a
% Bragg diffraction grating at an EUV wavelength. The grating comprises
% ruthenium lamellae supported on a thin, DLC (diamond-like carbon) film.
% (This configuration approximates the optical geometry at the periphery of
% a high-NA, diffractive EUV microlens.) The gdc_demo_engine parameters can
% be varied to represent a variety of slanted or unslanted lamellar grating
% structures.
%
% gdc_demo10 covers the following topics:
%   - parameterization
%   - coordinate break
%   - replication module
%   - refractive index tables
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
% >> gdc_demo10
% Elapsed time is 16.286061 seconds.
%  
% Wavelengths, diffraction efficiencies (unpolarized, order 2)
% 0.006   0.0047412
% 0.007   0.0013041
% 0.008    0.039929
% 0.009     0.22081
%  0.01     0.53283
% 0.011     0.62163
% 0.012     0.31418
% 0.013     0.04204
% 0.014     0.05458
% 0.015    0.051368
% 0.016    0.012174
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

demo=true; % true: run gdc_demo_engine; false: run gdc_engine

% Define parameters for grating structure and incident field. (Length
% quantities are in micron units.) The wavelength parameter is vectorized
% in dimension 1.
Bragg_wavelength=0.011; % EUV design wavelength for Bragg diffraction
Bragg_order=2; % diffraction order for Bragg diffraction
dev_angle=30*pi/180; % deviation angle for Bragg_order, Bragg_wavelength
wavelength=(0.006:0.0002:0.016).'; % wavelengths for efficiency calculation
t1=0.010; % base film thickness
t2=0.220; % grating layer thickness
wall_t=0.010; % grating wall thickness
% rep_count ("replication count") is the stratification number for the
% grating layer (preferably a power of 2).
rep_count=32;
m_max=20; % maximum diffraction order index

% Compute derived structure parameters. Functions of wavelength (film_pmt
% and grating_pmt) are size-matched to wavelength in dimension 1.
d=Bragg_order*Bragg_wavelength/sin(dev_angle); % grating period
slant=dev_angle/2; % slant angle of grating walls from the substrate normal
% Get tabulation wavelengths and refractive indices for the base film and
% for the grating layer; interpolate onto specified wavelengths.
[w,n]=read_nk('d-C.nk',0.0001); % DLC
n=interp1(w,n,wavelength);
film_pmt=n.^2; % base film permittivity
[w,n]=read_nk('Ru.nk',0.0001); % ruthenium
n=interp1(w,n,wavelength);
grating_pmt=n.^2; % grating layer permittivity
clear w n
sub_pmt=1; % substrate permittivity

if ~demo % The following code block is replicated in gdc_demo_engine.
    % Construct grating. (Note: Two grating periods must be specified
    % because the grating can generally be biperiodic, although for this
    % example the second period (d22,d32) is irrelevant.)
    clear grating
    grating.pmt={sub_pmt,film_pmt,grating_pmt,1}; % material permittivities
    grating.pmt_sub_index=1; % substrate permittivity index
    grating.pmt_sup_index=4; % superstrate permittivity index
    grating.d21=d; % first grating period: x2 projection
    grating.d31=0; % first grating period: x3 projection
    grating.d22=0; % second grating period: x2 projection
    grating.d32=d; % second grating period: x3 projection
    % Construct the stratum for the base film.
    clear stratum
    stratum.type=0; % homogeneous stratum
    stratum.pmt_index=2; % stratum's permittivity index
    stratum.thick=t1; % stratum thickness
    grating.stratum{1}=stratum;
    % Construct the stratum for the grating layer. First construct a cell
    % array of strata comprising the replication module.
    clear stratum
    stratum{1}.type=1; % module's first stratum is uniperiodic
    stratum{1}.thick=t2/rep_count; % module thickness
    % The following h11, h12 spec indicates that the stratum's period
    % vector matches the first grating period (GD-Calc.pdf, equations 3.22
    % and 3.23).
    stratum{1}.h11=1;
    stratum{1}.h12=0;
    clear stripe
    stripe.c1=0; % first stripe's boundary on positive side
    stripe.pmt_index=4; % first stripe's permittivity index
    stratum{1}.stripe{1}=stripe;
    stripe.c1=...
        wall_t/(d*cos(slant)); % second stripe's boundary on positive side
    stripe.pmt_index=3; % second stripe's permittivity index
    stratum{1}.stripe{2}=stripe;
    stratum{2}.type=3; % module's second stratum is a coordinate break.
    stratum{2}.dx2=-stratum{1}.thick*tan(slant); % x2 shift of upper strata
    stratum{2}.dx3=0; % x3 shift of upper strata
    % Construct the grating's second stratum as a replication module.
    grating.stratum{2}.type=4; % replication module
    grating.stratum{2}.stratum=stratum;
    grating.stratum{2}.rep_count=rep_count; % replication count
    clear stratum stripe
end

% Define the indicent field (normal incidence).
clear inc_field
inc_field.wavelength=wavelength;
inc_field.f2=0;
inc_field.f3=0;

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
        gdc_demo_engine(10,inc_field,order,...
        sub_pmt,film_pmt,grating_pmt,d,t1,t2,rep_count,wall_t,slant);
else
    [param_size,scat_field,inc_field]=gdc(grating,inc_field,order);
end
toc
if isempty(scat_field)
    disp('Interrupted by user.');
    return
end

% Compute the diffraction efficiencies. (Only the transmitted waves are
% retained.)
[R,T]=gdc_eff(scat_field,inc_field);
clear R
% Keep only the Bragg order.
T=T([T.m1]==Bragg_order);
% Tabulate and plot the first-order diffraction efficiency for unpolarized
% incident illumination.
disp(' ');
p1=1:5:length(wavelength); % Select wavelength tabulation indices.
disp(['Wavelengths, diffraction efficiencies (unpolarized, order ' ...
    num2str(Bragg_order) ')']);
disp(num2str([wavelength(p1), 0.5*([T.eff1(p1)]+[T.eff2(p1)])]));
figure, plot(wavelength, 0.5*([T.eff1]+[T.eff2]),'k');
title(['Efficiency in order ' num2str(Bragg_order)]);
xlabel('wavelength (micron)');
ylabel('eff');

% Plot the grating.
clear pmt_display
pmt_display(1).name='';
pmt_display(1).color=[];
pmt_display(1).alpha=1;
pmt_display(2).name='DLC';
pmt_display(2).color=[.75,.75,.75];
pmt_display(2).alpha=1;
pmt_display(3).name='Ru';
pmt_display(3).color=[.875,.875,.875];
pmt_display(3).alpha=1;
pmt_display(4).name='';
pmt_display(4).color=[];
pmt_display(4).alpha=1;
x_limit=[-0.5*(t1+t2),-(t1+t2),-(t1+t2);1.5*(t1+t2),(t1+t2),(t1+t2)];
gdc_plot(grating,1,pmt_display,x_limit);
