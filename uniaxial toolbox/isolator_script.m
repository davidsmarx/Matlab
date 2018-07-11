n1 = 1.0;
%no = 1.54426; ne = 1.55335; % calcite
no = 1.94657144; ne = 2.15093137; % YVO4
aoi = 10*pi/180; %2.3614*pi/180;
k1 = [0 sin(aoi) cos(aoi)]'; % for aoi in the x-y plane
ei = [1 1 0]'; ei = ei./sqrt(ei'*ei);

% first wedge:
c = [sin(22.5*pi/180) cos(22.5*pi/180) 0]'; % crystal axis
n = [0 sin(0*pi/180) cos(0*pi/180)]'; % surface normal
%
%[ke, se, me, ko] = uniaxial_raytracing_iso2uni(k1, n, c, n1, no, ne);
System/Prescription Data

File : H:\polarizing beam combiner\ray tracing\wollaston_test.ZMX
Title: Lens has no title.
Date : MON JAN 29 2001
Configuration 1 of 2

LENS NOTES:

       Notes...
        
        

GENERAL LENS DATA:

Surfaces                :               16
Stop                    :               11
System Aperture         : Float By Stop Size = 0.392013
Glass Catalogs          : Schott MY_GLASS_CATALOG
Ray Aiming              : Off
Apodization             :Uniform, factor =   0.00000E+000
Effective Focal Length  :       -1.421153 (in air)
Effective Focal Length  :       -1.421153 (in image space)
Back Focal Length       :        5.360254
Total Track             :        12.71292
Image Space F/#         :        2.126519
Paraxial Working F/#    :        2.477684
Working F/#             :         2.59338
Image Space NA          :       0.1978137
Object Space NA         :       0.1997414
Stop Radius             :      -0.3920135
Paraxial Image Height   :      0.07990391
Paraxial Magnification  :      -0.9987988
Entrance Pupil Diameter :       0.6683004
Entrance Pupil Position :         -2.3944
Exit Pupil Diameter     :        4.037762
Exit Pupil Position     :        10.00569
Field Type              : Object height in Millimeters
Maximum Field           :            0.08
Primary Wave            :           1.455
Lens Units              :   Millimeters
Angular Magnification   :       0.2383167

Fields          : 2
Field Type: Object height in Millimeters
#        X-Value        Y-Value         Weight
 1       0.000000       0.080000       1.000000
 2       0.000000      -0.080000       1.000000

Vignetting Factors
#       VDX       VDY       VCX       VCY
1  0.000000  0.000000  0.000000  0.000000
2  0.000000  0.000000  0.000000  0.000000

Wavelengths     : 3
Units: Microns
#          Value         Weight
 1       1.455000       1.000000
 2       1.420000       1.000000
 3       1.495000       1.000000

SURFACE DATA SUMMARY:

Surf     Type              Comment         Radius      Thickness                Glass      Diameter          Conic
 OBJ TILTSURF                                   -         1e-006               SILICA          0.16              -
   1 TILTSURF                                   -           0.24                          0.1600004              -
   2 GRINSUR9         GRIN @ 1.420       Infinity       4.432922              SLW-1.8           1.8              0
   3 STANDARD                            Infinity              1                                1.8              0
   4 COORDBRK    ROTATE WAVE PLATE              -              0                                  -              -
   5 JONESMAT      HALF-WAVE PLATE              -              1                          0.9844816              -
   6 COORDBRK          ROTATE BACK              -              0                                  -              -
   7 COORDBRK     VARIATION IN AOI              -              0                                  -              -
   8 BIRE__IN                            Infinity              1                 YVO4           1.5              0
   9 COORDBRK          WEDGE ANGLE              -              0                                  -              -
  10 BIRE_OUT                            Infinity              0                                1.5              0
 STO STANDARD                            Infinity            0.1                          0.7840269              0
  12 BIRE__IN                            Infinity              0                 YVO4           1.5              0
  13 COORDBRK          AFTER WEDGE              -              1                                  -              -
  14 BIRE_OUT                            Infinity              2                                1.5              0
  15 PARAXIAL                                   -           1.94                                  2              -
 IMA STANDARD                            Infinity                                         0.3792012              0

SURFACE DATA DETAIL:

Surface OBJ     : TILTSURF
 X Tangent      :                0
 Y Tangent      :           0.1405
 Scattering     : None
Surface   1     : TILTSURF
 X Tangent      :                0
 Y Tangent      :           0.1405
 Scattering     : None
Surface   2     : GRINSUR9
 Comment        : GRIN @ 1.420
 Delta t        :              0.1
 X Tangent      :                0
 Y Tangent      :           0.1405
 Aperture       : Floating Aperture
 Maximum Radius :           0.9
 Scattering     : None
Surface   3     : STANDARD
 Aperture       : Floating Aperture
 Maximum Radius :           0.9
 Scattering     : None
Surface   4     : COORDBRK
 Comment        : ROTATE WAVE PLATE
 Decenter X     :                0
 Decenter Y     :                0
 Tilt About X   :                0
 Tilt About Y   :                0
 Tilt About Z   :                0
 Order          : Decenter then tilt
 Scattering     : None
Surface   5     : JONESMAT
 Comment        : HALF-WAVE PLATE
 A real         :                1
 A imag         :                0
 B real         :                0
 B imag         :                0
 C real         :                0
 C imag         :                0
 D real         :               -1
 D imag         :                0
 Scattering     : None
Surface   6     : COORDBRK
 Comment        : ROTATE BACK
 Decenter X     :                0
 Decenter Y     :                0
 Tilt About X   :                0
 Tilt About Y   :                0
 Tilt About Z   :                0
 Order          : Tilt then decenter
 Scattering     : None
Surface   7     : COORDBRK
 Comment        : VARIATION IN AOI
 Decenter X     :                0
 Decenter Y     :                0
 Tilt About X   :                0
 Tilt About Y   :                0
 Tilt About Z   :                0
 Order          : Decenter then tilt
 Scattering     : None
Surface   8     : BIRE__IN
 Mode Flag      :                0
 X-Cosine       :                0
 Y-Cosine       :                1
 Z-Cosine       :                0
 Aperture       : Floating Aperture
 Maximum Radius :          0.75
 Scattering     : None
Surface   9     : COORDBRK
 Comment        : WEDGE ANGLE
 Decenter X     :                0
 Decenter Y     :                0
 Tilt About X   :        11.400241
 Tilt About Y   :                0
 Tilt About Z   :                0
 Order          : Decenter then tilt
 Scattering     : None
Surface  10     : BIRE_OUT
 Aperture       : Floating Aperture
 Maximum Radius :          0.75
 Scattering     : None
Surface STO     : STANDARD
 Scattering     : None
Surface  12     : BIRE__IN
 Mode Flag      :                1
 X-Cosine       :                1
 Y-Cosine       :                0
 Z-Cosine       :                0
 Aperture       : Floating Aperture
 Maximum Radius :          0.75
 Scattering     : None
Surface  13     : COORDBRK
 Comment        : AFTER WEDGE
 Decenter X     :                0
 Decenter Y     :                0
 Tilt About X   :       -11.400241
 Tilt About Y   :                0
 Tilt About Z   :                0
 Order          : Decenter then tilt
 Scattering     : None
Surface  14     : BIRE_OUT
 Aperture       : Floating Aperture
 Maximum Radius :          0.75
 Scattering     : None
Surface  15     : PARAXIAL
 Focal length   :             1.94
 OPD Mode       :                0
 Scattering     : None
Surface IMA     : STANDARD
 Scattering     : None

COATING DEFINITIONS:


EDGE THICKNESS DATA:

Surf         X-Edge         Y-Edge
 OBJ       0.000001       0.000001
   1       0.240000       0.355210
   2       4.432922       4.306472
   3       1.000000       1.000000
   4       0.000000       0.000000
   5       1.000000       1.000000
   6       0.000000       0.000000
   7       0.000000       0.000000
   8       1.000000       1.000000
   9       0.000000       0.000000
  10       0.000000       0.000000
 STO       0.100000       0.100000
  12       0.000000       0.000000
  13       1.000000       1.000000
  14       2.000000       2.000000
  15       1.940000       1.940000
 IMA       0.000000       0.000000

MULTI-CONFIGURATION DATA:

Configuration   1:

  1 Param 1     8 :             0 
  2 Param 1    12 :             1 
  3 Y-field     1 :          0.08 

Configuration   2:

  1 Param 1     8 :             1 
  2 Param 1    12 :             0 
  3 Y-field     1 :       -2.3614 

SOLVE AND VARIABLE DATA:

 Parameter  2 Surf   1  : Chief ray follow field 1, wavelength 0
 Semi Diameter   2      : Fixed
 Semi Diameter   3      : Fixed
 Parameter  5 Surf   6  : Pickup from 4 times -1.000000, plus 0.000000
 Semi Diameter   8      : Fixed
 Parameter  3 Surf   9  : Variable
 Semi Diameter  10      : Fixed
 Semi Diameter  11      : Fixed
 Semi Diameter  12      : Fixed
 Parameter  3 Surf  13  : Pickup from 9 times -1.000000, plus 0.000000
 Semi Diameter  14      : Fixed
 Semi Diameter  15      : Fixed

INDEX OF REFRACTION DATA:



Surf                  Glass    Temp    Pres      1.455000    1.420000    1.495000
   0                 SILICA   20.00    1.00    1.44514433  1.44554946  1.44467653
   1                          20.00    1.00    1.00000000  1.00000000  1.00000000
   2                SLW-1.8   20.00    1.00    1.59064502  1.59083690  1.59044202
   3                          20.00    1.00    1.00000000  1.00000000  1.00000000
   4              <CRD BRK>                    1.00000000  1.00000000  1.00000000
   5                          20.00    1.00    1.00000000  1.00000000  1.00000000
   6              <CRD BRK>                    1.00000000  1.00000000  1.00000000
   7              <CRD BRK>                    1.00000000  1.00000000  1.00000000
   8                   YVO4   20.00    1.00    1.94657144  1.94729303  1.94577675
              Extraordinary                    2.15093137  2.15185305  2.14992451
   9              <CRD BRK>                    1.94657144  1.94729303  1.94577675
  10                          20.00    1.00    1.00000000  1.00000000  1.00000000
  11                          20.00    1.00    1.00000000  1.00000000  1.00000000
  12                   YVO4   20.00    1.00    1.94657144  1.94729303  1.94577675
              Extraordinary                    2.15093137  2.15185305  2.14992451
  13              <CRD BRK>                    1.94657144  1.94729303  1.94577675
  14                          20.00    1.00    1.00000000  1.00000000  1.00000000
  15                          20.00    1.00    1.00000000  1.00000000  1.00000000
  16                          20.00    1.00    1.00000000  1.00000000  1.00000000

THERMAL COEFFICIENT OF EXPANSION DATA:

Surf                  Glass     TCE *10E-6
   0                 SILICA     0.00000000
   1                            0.00000000
   2                SLW-1.8     0.00000000
   3                            0.00000000
   4              <CRD BRK>     0.00000000
   5                            0.00000000
   6              <CRD BRK>     0.00000000
   7              <CRD BRK>     0.00000000
   8                   YVO4     4.43000000
   9              <CRD BRK>     4.43000000
  10                            0.00000000
  11                            0.00000000
  12                   YVO4     4.43000000
  13              <CRD BRK>     4.43000000
  14                            0.00000000
  15                            0.00000000
  16                            0.00000000

F/# DATA:

F/# calculations consider vignetting factors and ignore surface apertures.

             Wavelength:        1.455000            1.420000            1.495000    
 #                Field        Tan       Sag       Tan       Sag       Tan       Sag
 1    0.0000, 0.0800 mm:    2.6115    2.5923    2.6150    2.5957    2.6080    2.5887
 2   0.0000, -0.0800 mm:    2.6064    2.5948    2.6099    2.5983    2.6029    2.5912

GLOBAL VERTEX COORDINATES, ORIENTATIONS, AND ROTATION/OFFSET MATRICES:

Reference Surface: 1

Surf          R11           R12           R13             X
              R21           R22           R23             Y
              R31           R32           R33             Z

  0      1.000000      0.000000      0.000000      0.000000 
         0.000000      1.000000      0.000000      0.000000
         0.000000      0.000000      1.000000     -0.000001

  1      1.000000      0.000000      0.000000      0.000000 
         0.000000      1.000000      0.000000      0.000000
         0.000000      0.000000      1.000000      0.000000

  2      1.000000      0.000000      0.000000      0.000000 GRIN @ 1.420
         0.000000      1.000000      0.000000      0.000000
         0.000000      0.000000      1.000000      0.240000

  3      1.000000      0.000000      0.000000      0.000000 
         0.000000      1.000000      0.000000      0.000000
         0.000000      0.000000      1.000000      4.672922

  4      1.000000      0.000000      0.000000      0.000000 ROTATE WAVE PLATE
         0.000000      1.000000      0.000000      0.000000
         0.000000      0.000000      1.000000      5.672922

  5      1.000000      0.000000      0.000000      0.000000 HALF-WAVE PLATE
         0.000000      1.000000      0.000000      0.000000
         0.000000      0.000000      1.000000      5.672922

  6      1.000000      0.000000      0.000000      0.000000 ROTATE BACK
         0.000000      1.000000      0.000000      0.000000
         0.000000      0.000000      1.000000      6.672922

  7      1.000000      0.000000      0.000000      0.000000 VARIATION IN AOI
         0.000000      1.000000      0.000000      0.000000
         0.000000      0.000000      1.000000      6.672922

  8      1.000000      0.000000      0.000000      0.000000 
         0.000000      1.000000      0.000000      0.000000
         0.000000      0.000000      1.000000      6.672922

  9      1.000000      0.000000      0.000000      0.000000 WEDGE ANGLE
         0.000000      0.980270     -0.197661      0.000000
         0.000000      0.197661      0.980270      7.672922

 10      1.000000      0.000000      0.000000      0.000000 
         0.000000      0.980270     -0.197661      0.000000
         0.000000      0.197661      0.980270      7.672922

 11      1.000000      0.000000      0.000000      0.000000 
         0.000000      0.980270     -0.197661      0.000000
         0.000000      0.197661      0.980270      7.672922

 12      1.000000      0.000000      0.000000      0.000000 
         0.000000      0.980270     -0.197661     -0.019766
         0.000000      0.197661      0.980270      7.770949

 13      1.000000      0.000000      0.000000      0.000000 AFTER WEDGE
         0.000000      1.000000      0.000000     -0.019766
         0.000000      0.000000      1.000000      7.770949

 14      1.000000      0.000000      0.000000      0.000000 
         0.000000      1.000000      0.000000     -0.019766
         0.000000      0.000000      1.000000      8.770949

 15      1.000000      0.000000      0.000000      0.000000 
         0.000000      1.000000      0.000000     -0.019766
         0.000000      0.000000      1.000000     10.770949

 16      1.000000      0.000000      0.000000      0.000000 
         0.000000      1.000000      0.000000     -0.019766
         0.000000      0.000000      1.000000     12.710949


ELEMENT VOLUME DATA:

Values are only accurate for plane and spherical surfaces.
Element volumes are computed by assuming edges are squared up
to the larger of the front and back radial aperture.
Single elements that are duplicated in the Lens Data Editor
for ray tracing purposes may be listed more than once yielding
incorrect total mass estimates.

                             Volume cc   Density g/cc         Mass g
Element surf   2 to   3       0.011280       3.600000       0.040609
Total Mass:                                                 0.040609

CARDINAL POINTS:

Object space positions are measured with respect to surface 1.
Image space positions are measured with respect to the image surface.
The index in both the object space and image space is considered.

                           Object Space         Image Space
W = 1.455000 (Primary)
Focal Length      :            2.053550           -1.421000
Principal Planes  :           -4.108061            2.841253
Nodal Planes      :           -3.475511            3.473803
Focal Planes      :           -2.054512            1.420254
Anti-Nodal Planes :           -0.000962           -0.000746

W = 1.420000
Focal Length      :            2.052199           -1.419667
Principal Planes  :           -4.102288            2.839439
Nodal Planes      :           -3.469756            3.471971
Focal Planes      :           -2.050089            1.419772
Anti-Nodal Planes :            0.002110            0.000105

W = 1.495000
Focal Length      :            2.054910           -1.422401
Principal Planes  :           -4.114030            2.843156
Nodal Planes      :           -3.481522            3.475664
Focal Planes      :           -2.059120            1.420754
Anti-Nodal Planes :           -0.004211           -0.001647
[ke, se, me, ko, eetmp, eotmp] = uniaxial_raytracing_iso2uni_eig(k1, n, c, n1, no, ne);
[eo, ee, ers, erp, to, te, rs, rp] = uniaxial_match_iso2uni(k1,ko,ke,se,n,[n1;no;me],c,ei);
Eo = sqrt(no)*to*eo;
Ee = sqrt(ne)*te*ee;
Ers = sqrt(n1)*rs*ers;
Erp = sqrt(n1)*rp*erp;
return

% exit first wedge:
n = [0 0 1]';
% ext-beam:
kea = snells_law(ke,n,me,1.0); % ray vector in air
[kere, sere, mer, kero] = uniaxial_raytracing_uniref(ke,n,c,me,no,ne); % reflected e&o-rays
[ees, eep, eeo, eee, tes, tep, reo, ree] = ...
   uniaxial_match_uni2iso(ke,kea,kero,kere,sere,n,me,no,mer,1.0,c,ee);

% ord-beam:
koa = snells_law(ko,n,no,1.0);
[kore, sore, koro] = uniaxial_raytracing_uniref(ko,n,c,no,no,ne);

