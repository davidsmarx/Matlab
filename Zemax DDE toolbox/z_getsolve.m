function solvetype = z_getsolve(zchan, n_surf, code, pnum)
% solvetype = z_getsolve(zchan, n_surf, code, pnum)
% 
% This item returns data about any solve on any surface. The syntax is 
% GetSolve,surface,code
% 
% where code is an integer code indicating which surface parameter the solve
% data is for. The solve data is returned in the following
% formats, depending upon the code value.
%
% if code is parameter or extra data value, pnum indicates which parameter
% or extra data value
% 
% The solvetype is an integer code, and the parameters have meanings that
% depend upon the solve type; see
% the chapter Solves for details. See also SetSolve.
% 
% GetSolve Code Returned data format
% 0, curvature solvetype, parameter1, parameter2
% 1, thickness solvetype, parameter1, parameter2, parameter3
% 2, glass solvetype, pickupsurf
% 3, semi-diameter solvetype, pickupsurf
% 4, conic solvetype, pickupsurf
% 5-12, parameters 1-8 solvetype, pickupsurf, offset, scalefactor
% 1001+, extra data values
% 1+
% solvetype, pickupsurf, scalefactor

% validate code
if ischar(code)
    switch lower(code)
        case {'curv','curvature'}
            code = 0;
        case {'thick','thickness'}
            code = 1;
        case {'glass'}
            code = 2;
        case {'semi-diameter','sdia'}
            code = 3;
        case {'conic'}
            code = 4;
        case {'parameter','parm','param'},
            code = 4 + pnum;
        case {'extra data','extra'},
            code = 1000 + pnum;
        otherwise
            error(['code ' code ' is not recognized code type']);
    end
else
    if code < 0,
        error([num2str(code) ' is not a valid code']);
    end
end

% Chapter 13: SOLVES 349
% SUMMARY OF SOLVES
% Solve type First parameter Second parameter Third parameter Code Integer
% Curvature: Fixed 0
% Curvature: Variable V 1
% Curvature: Marginal ray
% angle
% Angle M 2
% Curvature: Chief ray angle Angle C 3
% Curvature: Pick up Surface Scale factor P 4
% Curvature: Marginal ray
% normal
% N 5
% Curvature: Chief ray normal N 6
% Curvature: Aplanatic A 7
% Curvature: Element power Power X 8
% Curvature: Concentric with
% surface
% Surface to be
% concentric to
% S 9
% Curvature: Concentric with
% radius
% Surface to be
% concentric with
% R 10
% Curvature: F/# Paraxial F/# F 11
% Thickness: Fixed 0
% Thickness: Variable V 1
% Thickness: Marginal ray
% height
% Height Pupil zone M 2
% Thickness: Chief ray height Height C 3
% Thickness: Edge thickness Thickness Radial height
% (use zero for
% semi-diameter)
% E 4
% Thickness: Pick up Surface Scale factor Offset P 5
% Thickness: Optical path
% difference
% OPD Pupil zone O 6
% Thickness: Position Surface Length f rom
% surface
% T 7
% Thickness: Compensator Surface Sum of surface
% thicknesses
% S 8
% Thickness: Center of
% curvature
% Surface to be at
% the COC of
% X 9
% Glass: Fixed 0
% Glass: Model Index Nd Abbe Vd Dpgf 1
% Glass: Pick up Surface P 2
% Glass: Substitute Catalog name S 3
% Glass: Offset Index Nd offset Abbe Vd offset O 4
% Semi-Diameter: Automatic 0
% Semi-Diameter: Fixed U 1
% Semi-Diameter: Pick up Surface P 2
% Semi-Diameter: Maximum M 3
% Conic: Fixed 0
% Conic: Variable V 1
% Conic: Pick up Surface Scale factor P 2
% Parameter: Fixed 0
% Parameter: Variable V 1
% Parameter: Pick up Surface Scale factor Offset P 2
% Parameter: Chief ray Field Wavelength C 3
% Solve type First parameter Second parameter Third parameter Code Integer

cmdstr = ['GetSolve,' num2str(n_surf) ',' num2str(code)]
solvetype = ddereq(zchan,cmdstr,[1 1]);

