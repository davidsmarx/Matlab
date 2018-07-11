cls
% 
% 

demo=false; % true: run gdc_demo_engine; false: run gdc_engine

clear grating
% substrate & superstrate = air
grating.pmt = {1, 3.695.^2};
grating.pmt_sub_index = 1;
grating.pmt_sup_index = 1;
% vector d1 = x-direction
d_g_1 = 1.5*[1; 0]; % period = 1.5
grating.d21 = d_g_1(1); grating.d31 = d_g_1(2); 
% vector d2 = y-direction
d_g_2 = 5*[0; 1]; % period = 5
grating.d22 = d_g_2(1); grating.d32 = d_g_2(2);

% two strata: uniperiodic grating, homogeneous wafer layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stratum #1: uniperiodic grating, 50% duty cycle
clear stratum;
stratum.type = 1;
stratum.thick = 25;
% grating period and direction
d_s = d_g_1;
hh = [d_g_1./(d_g_1'*d_g_1) d_g_2./(d_g_2'*d_g_2)] \ (d_s./(d_s'*d_s));
stratum.h11 = hh(1); stratum.h12 = hh(2);
% this stratum has two stripes
stratum.stripe{1} = struct('c1', 0.25, 'pmt_index',2);
stratum.stripe{2} = struct('c1', 0.75, 'pmt_index',1);
grating.stratum{1} = stratum;

% stratum #2: homogeneous layer of silicon
clear stratum;
stratum.type = 0;
stratum.thick = 600;
stratum.pmt_index = 2;
grating.stratum{2} = stratum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% incident field
% Define the indicent field.
clear inc_field
inc_field.wavelength=1.31;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector. The grating-normal (x1)
% projection is implicitly f1=-cos(phi)/wavelength.
theta = 5*P; phi = 10*P;
inc_field.f2=sin(theta)*cos(phi)/inc_field.wavelength;
inc_field.f3=sin(theta)*sin(phi)/inc_field.wavelength;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diffraction orders
% for each diffraction order along vector d2, list the orders along d1
order.m2 = 0;
order.m1 = -10:10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the call
gdc(grating,inc_field,order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the grating
% pmt_display => one for each entry in grating.pmt{}
%grating.pmt = {1, 3.695.^2};
pmt_display(1) = struct('name','','color',[],'alpha',1);
pmt_display(2) = struct('name','Si','color',[0 0 1],'alpha',1);
x_limit=[
    -700 -4 -4
    700 4 4
    ];
h_plot=gdc_plot(grating,1,pmt_display,x_limit);




return

% Define parameters for grating structure:
d1=1.50; % grating period in x2 direction
d2=1.00; % grating period in x3 direction
h=0.25; % grating height
grating_pmt=1.5^2; % grating permittivity
L1=16; % number of grating strata
m_max=5; % maximum diffraction order index

% incident field
wavelength=1.533;
phi=30*pi/180; % incidence polar angle from x1 axis
psi=45*pi/180; % incidence azimuthal angle around x1 axis (zero on x2 axis)

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
    %[param_size,scat_field,inc_field]=gdc(grating,inc_field,order);
    gdc(grating,inc_field,order)
end
toc

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

return

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

