% Code examples from GD-Calc_Intro.pdf  (software.kjinnovation.com)
% Version 02/13/2006
%
% These examples require that gdc.m and gdc_plot.m be on the MATLAB path.
% (To run the examples interactively, select the code, right click, and
% choose "Evaluate Selection".)

% Listing (6)
d=1.5; % grating period
grating_pmt=(1.52+6.46*i)^2; % grating permittivity
clear grating
grating.pmt={1.0,grating_pmt};
grating.pmt_sub_index=2;
grating.pmt_sup_index=1;
grating.d21=d;
grating.d31=0;
grating.d22=0;
grating.d32=d;
grating.stratum={};

% Data validation
gdc(grating);

% Listing (9)
thick=0.5; % stratum thickness
grating.pmt_sub_index=1;
clear stratum
stratum.type=0; % homogeneous
stratum.thick=thick;
stratum.pmt_index=2;
grating.stratum{1}=stratum;

% Data validation
gdc(grating);

% Listing (10)
width=0.5; % rod width
clear stratum
stratum.type=1; % uniperiodic
stratum.thick=thick;
stratum.h11=1;
stratum.h12=0;
clear stripe
stripe.c1=-0.5*width/d;
stripe.pmt_index=1;
stratum.stripe{1}=stripe;
stripe.c1=0.5*width/d;
stripe.pmt_index=2;
stratum.stripe{2}=stripe;
grating.stratum{1}=stratum;

% Listing (11)
clear pmt_display
pmt_display(1).name='';
pmt_display(1).color=[];
pmt_display(1).alpha=1;
pmt_display(2).name='Tungsten';
pmt_display(2).color=[1,1,1]*0.75;
pmt_display(2).alpha=1;
x_limit=[-thick,-1.75*d,-1.75*d;...
    2*thick,1.75*d,1.75*d];
gdc_plot(grating,1,pmt_display,x_limit);

% Listing (12)
clear stratum
stratum.type=2; % biperiodic
stratum.thick=thick;
stratum.h11=1;
stratum.h12=0;
stratum.h21=0;
stratum.h22=1;
clear stripe
stripe.type=1; % inhomogeneous
stripe.c1=-0.5*width/d;
clear block
block.c2=-0.5*width/d;
block.pmt_index=1;
stripe.block{1}=block;
block.c2=0.5*width/d;
block.pmt_index=2;
stripe.block{2}=block;
stratum.stripe{1}=stripe;
clear stripe
stripe.type=0; % homogeneous
stripe.c1=0.5*width/d;
stripe.pmt_index=2;
stratum.stripe{2}=stripe;
grating.stratum{1}=stratum;

% Plot the grating.
gdc_plot(grating,1,pmt_display,x_limit);

% Listing (18)
clear stratum
stratum.type=1; % uniperiodic
stratum.thick=thick;
stratum.h11=1;
stratum.h12=0;
clear stripe
stripe.c1=-0.5*width/d;
stripe.pmt_index=1;
stratum.stripe{1}=stripe;
stripe.c1=0.5*width/d;
stripe.pmt_index=2;
stratum.stripe{2}=stripe;
grating.stratum{1}=stratum; % same as listing (10)
stratum.h11=0; % Swap h11 and h12.
stratum.h12=1;
grating.stratum{2}=stratum;

% Plot the grating.
x_limit(2,1)=3*thick;
gdc_plot(grating,1,pmt_display,x_limit);

% Listing (19)
clear stratum
stratum.type=3; % coordinate break
stratum.dx2=d/2; % half-period x2-shift
stratum.dx3=d/2; % half-period x3-shift
grating.stratum{3}=stratum;
grating.stratum{4}=grating.stratum{1};
grating.stratum{5}=grating.stratum{2};

% Plot the grating.
x_limit(2,1)=5*thick;
gdc_plot(grating,1,pmt_display,x_limit);

% Listing (20)
for l1=4:11
    grating.stratum{l1}=grating.stratum{l1-3};
end

% Plot the grating.
x_limit(2,1)=9*thick;
gdc_plot(grating,1,pmt_display,x_limit);

% Listing (21)
clear stratum
stratum.type=4; % replication module
stratum.stratum{1}=grating.stratum{1};
stratum.stratum{2}=grating.stratum{2};
stratum.stratum{3}=grating.stratum{3};
stratum.rep_count=4; % replication count
grating.stratum={stratum};

% Plot the grating.
gdc_plot(grating,1,pmt_display,x_limit);

% Listing 24
m_max=10;
order=[];
m1=-m_max:m_max;
for m2=-m_max:m_max
    order(end+1).m2=m2;
    order(end).m1=m1;
end

% Display the order struct.
order.m2
order.m1

% Listing 25
m_max=10;
order=[];
m1=-m_max:m_max;
for m2=-m_max:m_max
    order(end+1).m2=m2;
    order(end).m1=m1(abs(m1)+abs(m2)<=m_max);
end

% Display the order struct.
order.m2
order.m1

% Listing 26
m_max=10;
order=[];
m1=-m_max:m_max;
for m2=-m_max:m_max
    order(end+1).m2=m2;
    order(end).m1 = m1(mod(m1-m2,2)==0);
end

% Display the order struct.
order.m2
order.m1
