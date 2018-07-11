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

set(gca,'dataaspectRatio',[0.01 1 1])

