% mks system
global CM MM UM NM IN MIL FS PS NS US MS MIN HH DAY P PARCMIN PARCSEC;
global FAHR2K FAHR2C CELS2K;

CM = 1e-2; % centimeters
MM = 1e-3; % millimeters
UM = 1e-6; % micrometers
NM = 1e-9; % nanometers
PM = 1e-12; % picometers
IN = 0.0254; % inches
MIL = 1e-3*IN; % mils
UIN = 1e-6*IN;

FS = 1e-15; % femto seconds
PS = 1e-12; %
NS = 1e-9;
US = 1e-6; % microseonds
MS = 1e-3; % milliseconds
MIN = 60;  % minutes
HH  = 60*60; % hours
DAY = 60*60*24; % days

HZ  = 1;
KHZ = 1e3; % kilohertz
MHZ = 1e6; % megahertz
GHZ = 1e9;
THZ = 1e12;
WAVENUMBER = 100*2.99792458e8; % wavenumber in 1/cm

V = 1;
MV = 1e-3; % mV

P = pi/180; % degrees
PARCMIN = P/60;
PARCSEC = PARCMIN/60;
MRAD = 1e-3;

% pressure
global PASCAL GPA DYN_CM2 PSI;
PASCAL = 1; % N/m^2
GPA = 1e9;
DYN_CM2 = 0.1; % dyn/cm^2
PSI = 4.448222/(IN*IN); % 1 psi = 6894.757 Pa

% Temperature
FAHR2K = @(f)((5./9).*(f-32) + 273.15);
FAHR2C = @(f)((5./9).*(f-32));
CELS2K = @(c)(c + 273.15);
C2FAHR = @(c)((9./5)*c + 32);

% Star Magnitude
STARMAG = 100.^(1/5);
