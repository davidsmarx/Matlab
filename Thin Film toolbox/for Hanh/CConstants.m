classdef CConstants
    
    properties (Constant = true)
        %%%%%%%%%%%%%%%%%%%% MKS UNITS %%%%%%%%%%%%%%%%%%%%%%%%%
        CM = 1e-2; % centimeters
        MM = 1e-3; % millimeters
        UM = 1e-6; % micrometers
        NM = 1e-9; % nanometers
        PM = 1e-12; % picometers
        IN = 0.0254; % inches
        MIL = 1e-3*0.0254; % mils
        UIN = 1e-6*0.0254;
        
        FS = 1e-15; % femto seconds
        PS = 1e-12; %
        NS = 1e-9;
        US = 1e-6; % microseonds
        MS = 1e-3; % milliseconds
        MIN = 60;  % minutes
        HH  = 60*60; % hours
        DAY = 60*60*24; % days
        
        MHZ = 1e6; % megahertz
        GHZ = 1e9;
        THZ = 1e12;
        WAVENUMBER = 100*2.99792458e8; % wavenumber in 1/cm
        
        P = pi/180; % degrees
        DEG = pi/180;
        PARCMIN =  pi/180/60;
        PARCSEC = pi/180/60/60;

        % pressure
        PASCAL = 1; % N/m^2
        GPA = 1e9;
        DYN_CM2 = 0.1; % dyn/cm^2
        PSI = 4.448222/( 0.0254* 0.0254); % 1 psi = 6894.757 Pa

        %%%%%%%%%%%%%%%%%%%% PHYSICAL CONSTANTS %%%%%%%%%%%%%%%%%%
        C = 2.99792458e8; % [vacuum m/s]
        E = exp(1);
        MU0 = 4*pi*1e-7; % [vacuum permeability N / A^2]
        E0 = 8.85418781e-012 % =1./(MU0*C*C); % [vacuum permittivity ]
        ECHARGE = 1.60217653e-19; % [electron charge Coulombs]
        HPLANCK = 6.6260693e-34; % [J s]
        KBOLTZMANN = 1.3806505e-23; % [J / K]
        g_n = 9.80665; % [m / s^2], from http://www.nist.gov/calibrations/vibration_measurements.cfm
    end % properties
    
    methods
        % Temperature
        function k = FAHR2K(S,f)
            k = ((5./9).*(f-32) + 273.15);
        end
        function c = FAHR2C(S,f)
            c = ((5./9).*(f-32));
        end
        function k = CELS2K(S,c)
            k = c + 273.15;
        end
        function f = C2FAHR(S,c)
            f = ((9./5)*c + 32);
        end
    end
    
end % classdef