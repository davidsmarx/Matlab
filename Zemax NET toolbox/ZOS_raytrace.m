classdef ZOS_raytrace

    properties

        Craytrace

    end % properties

    properties (Constant)
        MM = 1e-3;

    end % constants

    methods

        function S = ZOS_raytrace(Zapp)
            
            S.Craytrace = Zapp.PrimarySystem.Tools.OpenBatchRayTrace;

        end % init function

        function [vignetteCode, raypos, raycos, surfnorm, intensity, opd] = SingleRayNormUnpol(S, surfnumber, wavenumber, h, p)
            % [vignetteCode, raypos, raycos, surfnorm, intensity, opd] = SingleRayNormUnpol(surfnumber, wavenumber, h, p)
            % replaces z_gettrace
            % h = [hx hy]; p = [px py];
            
            % [ErrorCode, vignetteCode, xo, yo, zo, lo, mo, no, l2o, m2o, n2o, opd, intensity] = raytrace.SingleRayNormUnpol(ZOSAPI.Tools.RayTrace.RaysType.Real, Lens.Lens.OPTICAL_STOP.SurfaceNumber, 0, 0, 0, 0, 0, false);
            [retbool, ErrorCode, vignetteCode, xo, yo, zo, lo, mo, no, l2o, m2o, n2o, opd, intensity] = ...
                S.Craytrace.SingleRayNormUnpol(ZOSAPI.Tools.RayTrace.RaysType.Real, surfnumber, wavenumber, h(1), h(2), p(1), p(2), false);

            %disp(retbool);

            % evaluate error code here rather than return it
            if ErrorCode ~= 0,
                warning(['SingleRayNormUnpol Error Code: ' num2str(ErrorCode)]);
            end

            % return values
            raypos = [xo yo zo]'*S.MM;
            raycos = [lo mo no]';
            surfnorm = [l2o m2o n2o]';

        end % SingRayNormUnpol

        function [vignetteCode, raypos, raycos, surfnorm, rayI] = SingleRayDirectUnpol(S, wavenumber, nsurf_start, nsurf_stop, raypos_start, raycos_start)
            % [vignetteCode, raypos, raycos, surfnorm, rayI] = SingleRayDirectUnpol(S, wavenumber, nsurf_start, nsurf_stop, raypos_start, raycos_start)
            % replaces z_gettracedirect
            %
            % raypos_start in (m)

            % Zemax is in mm
            raypos_start = raypos_start/S.MM;

            % [ErrorCode, vignetteCode, xo, yo, zo, lo, mo, no, l2o, m2o, n2o, opd, intensity] = raytrace.SingleRayNormUnpol(ZOSAPI.Tools.RayTrace.RaysType.Real, Lens.Lens.OPTICAL_STOP.SurfaceNumber, 0, 0, 0, 0, 0, false);
            [retbool, ErrorCode, vignetteCode, xo, yo, zo, lo, mo, no, l2o, m2o, n2o, rayI] = ...
                S.Craytrace.SingleRayDirectUnpol(ZOSAPI.Tools.RayTrace.RaysType.Real, nsurf_start, nsurf_stop, wavenumber, ...
                raypos_start(1), raypos_start(2), raypos_start(3), raycos_start(1), raycos_start(2), raycos_start(3));

            %disp(retbool);

            % evaluate error code here rather than return it
            if ErrorCode ~= 0,
                warning(['SingleRayDirectUnpol Error Code: ' num2str(ErrorCode)]);
            end

            % return values
            raypos = [xo yo zo]'*S.MM;
            raycos = [lo mo no]';
            surfnorm = [l2o m2o n2o]';

        end % SingRayNormUnpol

    end % methods


end % classdef
    
    %function raytrace = ZOS_raytrace(Zapp)
% for now, just a how to example

% replacement for DDE GetTraceDirect
%raytrace = Zapp.PrimarySystem.Tools.OpenBatchRayTrace;

% When opening the tool, the user selects the number of maximum rays to be traced, the type of rays (the ray tracing can use real or paraxial rays), and the first and last surfaces where the rays will be traced.
% CreateDirectUnpol (int MaxRays, RaysType rayType, int startSurface, int toSurface)
% 
% The rays are then defined one by one in the AddRay function. The rays are defined by the wavelength and x, y, z, l, m, and n coordinates on any starting surface.
% AddRay (int waveNumber, double X, double Y, double Z, double L, double M, double N)
% 
% Results
% 
% This method can return:
% ReadNextResult (out int rayNumber, out int ErrorCode, out int vignetteCode, out double X, out double Y, out double Z, out double L, out double M, out double N, out double l2, out double m2, out double n2, out double intensity)
% •rayNumber
% •ErrorCode
% •vignetteCode: indicates if the ray was vignetted
% •X, Y, Z: the coordinates of the ray at the requested surface
% •L, M, N: direction cosines of the ray at the requested surface
% •l2, m2, n2: vector normal of the ray at the requested surface.
% •Intensity: the power of the ray (see discussion on the polarization in the IRayTraceNormUnpolData)
% 
% Next:


% prototype: CreateNormUnpol(ZemaxUI.ZOSAPI.Tools.BatchRayTraceTool this, int32 scalar maxRays, ZOSAPI.Tools.RayTrace.RaysType rayType, int32 scalar toSurface)
% raytracenorm = raytrace.CreateNormUnpol(5, ZOSAPI.Tools.RayTrace.RaysType.Real, dm1Surface);
% 
% startSurface = 27;
% toSurface = 55;
% raytracedirect = raytrace.CreateDirectUnpol(5, ZOSAPI.Tools.RayTrace.RaysType.Real, startSurface, toSurface);
% list_
% 
% % raypos = (x,y,z) (m) (assumes ZEMAX lens units = mm)
% % raycos = ray angles (cosx, cosy, cosz)
% raytracedirect.AddRay(1, 0, 0, 0, 0.001, 0, 0);

%% single ray trace
% SingleRayNormUnpol(), also SingleRayDirectUpnpl()
% 
% 
% 
% bool ZOSAPI.Tools.RayTrace.IBatchRayTrace.SingleRayNormUnpol  ( 
%   RaysType  rayType,  (e.g. ZOSAPI.Tools.RayTrace.RaysType.Real)
%   int  toSurf,  
%   int  waveNumber,  
%   double  Hx,  
%   double  Hy,  
%   double  Px,  
%   double  Py,  
%   bool  calcOPD,  
%   out int  ErrorCode,  
%   out int  vignetteCode,  
%   out double  xo,  
%   out double  yo,  
%   out double  zo,  
%   out double  lo,  
%   out double  mo,  
%   out double  no,  
%   out double  l2o,  
%   out double  m2o,  
%   out double  n2o,  
%   out double  opd,  
%   out double  intensity  
%  ) 
% [ErrorCode, vignetteCode, xo, yo, zo, lo, mo, no, l2o, m2o, n2o, opd, intensity] = raytrace.SingleRayNormUnpol(ZOSAPI.Tools.RayTrace.RaysType.Real, Lens.Lens.OPTICAL_STOP.SurfaceNumber, 0, 0, 0, 0, 0, false);
