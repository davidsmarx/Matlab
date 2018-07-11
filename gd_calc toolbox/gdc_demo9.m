% gdc_demo9: alignment sensor
%
% The grating structure in this demo comprises a phase-plate transmission
% grating in close proximity to a reflection grating. (Both elements are
% uniperiodic, lamellar gratings.) The lateral dislacement of the phase
% plate and the air gap between the gratings are both variable parameters,
% and the ratio of the +1 to -1 diffraction order efficiencies provides a
% sensitive measure of the lateral displacement.
%
% gdc_demo9 covers the following topics:
%	- multilayer grating (uniperiodic and homogeneous layers)
%   - harmonic indices
%   - parameterization
%   - coordinate break
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
% >> gdc_demo9
% Elapsed time is 4.980209 seconds.
%  
% Diffraction efficiencies for dx2=d*[0 0.125 0.25 0.375]
% R for t2=3*wavelength:
% -3    0.031538    0.047847    0.040499    0.002933
% -2   0.0086724   0.0088604    0.010334   0.0062346
% -1     0.11308     0.01542     0.10177     0.23084
%  0    0.018459    0.012857    0.020868    0.012857
%  1     0.11308     0.23084     0.10177     0.01542
%  2   0.0086724   0.0062346    0.010334   0.0088604
%  3    0.031538    0.002933    0.040499    0.047847
% R for t2=2*wavelength:
% -3    0.039677    0.047406    0.051274   0.0050904
% -2   0.0080337   0.0056262   0.0068248   0.0053821
% -1    0.097913    0.011035    0.093417     0.22771
%  0     0.01517    0.012441    0.017748    0.012441
%  1    0.097913     0.22771    0.093417    0.011035
%  2   0.0080337   0.0053821   0.0068248   0.0056262
%  3    0.039677   0.0050904    0.051274    0.047406
% R for t2=1*wavelength:
% -3    0.043042    0.037406    0.064805   0.0049872
% -2   0.0086529   0.0055333   0.0047146   0.0044991
% -1    0.092943   0.0056958    0.087841     0.25207
%  0    0.014058   0.0096163     0.01383   0.0096163
%  1    0.092943     0.25207    0.087841   0.0056958
%  2   0.0086529   0.0044991   0.0047146   0.0055333
%  3    0.043042   0.0049872    0.064805    0.037406
%
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

demo=true; % true: run gdc_demo_engine; false: run gdc_engine

% Define parameters for grating structure and incident field. The air space
% parameter (t2) is vectorized in the first dimension, and the phase
% plate's lateral shift parameter (dx2) is vectorized in the second
% dimension.
superstrate_pmt=2.25; % superstrate (phase plate) permittivity
substrate_pmt=16.0; % substrate permittivity
wavelength=0.633;
t1=wavelength/4; % bottom (reflecting) grating thickness
t2=wavelength*[3;2;1]; % air space thickness
t3=(wavelength/8)/(sqrt(superstrate_pmt)-1); % phase plate thickness
d=10*wavelength; % grating period
dx2=d*(0:63)/64; % phase plate's lateral shift
m_max=20; % maximum diffraction order index

if ~demo % The following code block is replicated in gdc_demo_engine.
    % Construct grating. (Note: Two grating periods must be specified
    % because the grating can generally be biperiodic, although for this
    % example the second period (d22,d32) is irrelevant.)
    clear grating
    grating.pmt=...
        {1,substrate_pmt,superstrate_pmt}; % material permittivities
    grating.pmt_sub_index=2; % substrate permittivity index
    grating.pmt_sup_index=3; % superstrate permittivity index
    grating.d21=d; % first grating period: x2 projection
    grating.d31=0; % first grating period: x3 projection
    grating.d22=0; % second grating period: x2 projection
    grating.d32=d; % second grating period: x3 projection
    % Construct the stratum for the reflecting grating.
    clear stratum
    stratum.type=1; % uniperiodic stratum
    stratum.thick=t1; % stratum thickness
    % The following h11, h12 spec indicates that the stratum's period
    % vector matches the first grating period (GD-Calc.pdf, equations 3.22
    % and 3.23).
    stratum.h11=1;
    stratum.h12=0;
    clear stripe
    stripe.c1=-0.25; % first stripe's boundary on positive side
    stripe.pmt_index=1; % first stripe's permittivity index
    stratum.stripe{1}=stripe;
    stripe.c1=0.25; % second stripe's boundary on positive side
    stripe.pmt_index=2; % second stripe's permittivity index
    stratum.stripe{2}=stripe;
    grating.stratum{1}=stratum;
    % Construct the stratum for the air space.
    clear stratum
    stratum.type=0; % homogeneous stratum
    stratum.pmt_index=1; % stratum's permittivity index
    stratum.thick=t2; % stratum thickness
    grating.stratum{2}=stratum;
    % Construct a (zero-thickness) stratum representing a coordinate break
    % (for the phase plate's lateral shift)
    clear stratum
    stratum.type=3; % coordinate break
    stratum.dx2=dx2; % x2 shift of all strata above the coordinate break
    stratum.dx3=0; % x3 shift of all strata above the coordinate break
    grating.stratum{3}=stratum;
    % Construct the stratum for the phase-plate grating.
    clear stratum
    stratum.type=1; % uniperiodic stratum
    stratum.thick=t3; % stratum thickness
    % The following h11, h12 spec indicates that the stratum's period
    % vector is half the first grating period (i.e., its fundamental
    % spatial frequency is twice that of the grating; see GD-Calc.pdf,
    % equations 3.22 and 3.23).
    stratum.h11=2;
    stratum.h12=0;
    clear stripe
    stripe.c1=-0.25; % first stripe's boundary on positive side
    stripe.pmt_index=1; % first stripe's permittivity index
    stratum.stripe{1}=stripe;
    stripe.c1=0.25; % second stripe's boundary on positive side
    stripe.pmt_index=3; % second stripe's permittivity index
    stratum.stripe{2}=stripe;
    grating.stratum{4}=stratum;
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
        gdc_demo_engine(9,inc_field,order,...
        substrate_pmt,superstrate_pmt,d,t1,t2,t3,dx2);
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
% grating (i.e., evanescent waves).
R=R(imag([scat_field.f1r])==0);
% Keep only orders -3...3.
R=R(abs([R.m1])<=3);
% Extract the diffraction order indices for the reflected waves.
R_m1=[R.m1].';
% Construct strings for plot legend.
m1_str={};
for k=1:length(R_m1)
    m1_str{k}=['m_1=' num2str(R_m1(k))];
end
% Tabulate and plot the diffraction order indices and diffraction
% efficiencies for unpolarized incident illumination.
disp(' ');
p2=1:ceil(length(dx2)/8):(length(dx2)/2); % Select dx2 tabulation indices.
disp(['Diffraction efficiencies for dx2=d*' mat2str(dx2(p2)/d)]);
for p1=1:length(t2)
    disp(['R for t2=' num2str(t2(p1)/wavelength) '*wavelength:']);
    R_eff=[];
    for k=1:length(R_m1)
        R_eff(k,:)=(R(k).eff1(p1,:)+R(k).eff2(p1,:))/2;
    end
    disp(num2str([R_m1 R_eff(:,p2)]));
    figure, plot(dx2.'/d,R_eff.');
    title(['Efficiency vs dx2 for t2=' ...
        num2str(t2(p1)/wavelength) '*wavelength']);
    % Note: "m1_str{:}" is MATLAB syntax for "m1_str{1},m1_str{2},...".
    legend(m1_str{:});
    xlabel('dx_2/d');
    ylabel('eff');
    axis([0,1,0,0.3]);
end

% Generate an animation showing a dx2 parameter scan.
clear pmt_display
pmt_display(1).name='';
pmt_display(1).color=[];
pmt_display(1).alpha=1;
pmt_display(2).name='reflector';
pmt_display(2).color=[1,1,1]*0.75;
pmt_display(2).alpha=1;
pmt_display(3).name='phase plate';
pmt_display(3).color=[1,1,1]*0.875;
pmt_display(3).alpha=1;
x_limit=[-0.5*(t1+t2(1)+t3),-1.5*d,-1.5*d;1.5*(t1+t2(1)+t3),1.5*d,1.5*d];
disp('Making movie, please wait ...');
F=struct('cdata',{},'colormap',{});
figure
p1=1;
for p2=1:4:length(dx2)
    % Plot the grating with air space t2(p1) and lateral shift dx2(p2).
    % (The "false" argument makes the plot go to the current figure.)
    gdc_plot(grating,[p1,p2],pmt_display,x_limit,false);
    F(end+1)=getframe;
end
disp('Done. Press any key to play movie.');
pause
movie(F,5);
