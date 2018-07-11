function [beamout, curvout, dxout, dyout] = ...
    FresnelPropagate(beam, dx, dy, lam, z, refindex, dxo, dyo, applycurv)
% [beamout, curvout, dxout, dyout] = ...
%    FresnelPropagate(beam, dx, dy, lam, z, refindex, dxo, dyo,
%    applycurv)
%
% modified Fresnel section of propagate()


% validate inputs
if abs(z) < lam, beamout = beam; return, end
if ~exist('dxo','var'), dxo = []; end
if ~exist('dyo','var'), dyo = []; end

%scalefactor = 1 + z*curv;
scalefactor = 1;

[Ny, Nx] = size(beam);

Gx = Nx*dx; Gy = Ny*dy;

xin = (-Nx/2:Nx/2-1)'*dx;
yin = (-Ny/2:Ny/2-1)'*dy;

% calculate Fresnel numbers based on grid size and divergence or
% convergence of initial wavefront
lambda_z = z * lam;
Fren_x = abs( scalefactor * Gx.^2 / lambda_z );
Fren_y = abs( scalefactor * Gy.^2 / lambda_z );

if ( Fren_x <= Nx & Fren_y <= Ny ),
    % small Fresnel number case => far field
    
    error('small Fresnel number, use propagate()');
    
elseif ( Fren_x > Nx & Fren_y > Ny ),
    % large Fresnel number case => near field
    
    % check scale of output grid to input grid, and apply geometric
    % curvature as necessary
    if ~isempty(dxo) & ~isempty(dyo),
        if abs( dxo/dx - dyo/dy ) <= 1e-6, % same scale
            scalefactor = sign(scalefactor) * dxo/dx;
                        
            geomCurv    = ( scalefactor - 1 ) ./ z;
            %beam = ApplyCurvature(beam, xin, yin, curv - geomCurv, lam);
            beam = ApplyCurvature(beam, xin, yin, -geomCurv, lam);
            dxout = dxo;
            dyout = dyo;
        else,
            error('Propagate: grid scale factor error');
        end
    else,
        error('not implemented');
        %geomCurv = curv;
        %dxout = abs(scalefactor) * dx;
        %dyout = abs(scalefactor) * dy;
    end
    
    lambda_z_eff = lambda_z./scalefactor;
    
    beamout = Fresnel; 
    
    curvout = geomCurv./scalefactor;
    xout = (-Nx/2:Nx/2-1)'*dxout;
    yout = (-Ny/2:Ny/2-1)'*dyout;

else,
    error('Propagate: Fresnel number conflict');
end
    
if exist('applycurv','var') & applycurv,
    beamout = ApplyCurvature(beamout, xout, yout, curvout, lam);
    curvout = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function beamout = Fresnel %(beam,dx,dy,lambda_z_eff,scalefactor)

        dfx = 1./Gx;
        dfy = 1./Gy;
        fx  = [-Nx/2:Nx/2-1]'*dfx;
        fy  = [-Ny/2:Ny/2-1]'*dfy;
        
        % calculate plane-wave spectrum of wavefront
        beam = fftshift(beam);
        beam = fft2(beam);
        beam = fftshift(beam);

        propconst = -i*pi*lambda_z_eff;

        efx = (dx*dfx./scalefactor)*exp( propconst * fx.^2 );
        efy = (dy*dfy)             *exp( propconst * fy.^2 );
        beam = (efy * efx.') .* beam;

        beam = fftshift(beam);
        if scalefactor < 0, 
            % propagating through focus inverts the wavefront
            beam = fft2(beam);
        else,
            %beam = ifft2(beam);
            beam = conj(fft2(conj(beam)));
        end
        beamout = fftshift(beam);

    end % Fresnel

end % main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beamout = ApplyCurvature(beam,x,y,curv,lambda)
% applies a curvature (quadratic phase) to a beam
% x and y are coordinate vectors

curvcoeff = pi*curv./lambda;

phase = exp(j*curvcoeff*y(:).^2) * exp(j*curvcoeff*x(:)'.^2);
beamout = phase .* beam;

clear phase;

end % ApplyCurvature
