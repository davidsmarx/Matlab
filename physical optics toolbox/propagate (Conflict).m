function [beamout, curvout, dxout, dyout] =...
    propagate(beam, dx, dy, curv, lam, z, dxo, dyo, applycurv)
% [beamout, curvout, dxout, dyout] =...
%    propagate(beam, dx, dy, curv, lam, z)
% [beamout, curvout, dxout, dyout] =...
%    propagate(beam, dx, dy, curv, lam, z, dxo, dyo, applycurv)
%
% replicates JPL FORTRAN code

% validate inputs
if abs(z) < lam, beamout = beam; curvout = curv; dxout = dxo; dyout = dyo; return, end
if ~exist('dxo','var'), dxo = []; end
if ~exist('dyo','var'), dyo = []; end
if ~exist('applycurv','var'), applycurv = true; end

scalefactor = 1 + z*curv;

[Ny, Nx] = size(beam);

Gx = Nx*dx; Gy = Ny*dy;

xin = (-Nx/2:Nx/2-1)'*dx;
yin = (-Ny/2:Ny/2-1)'*dy;

% calculate Fresnel numbers based on grid size and divergence or
% convergence of initial wavefront
lambda_z = z * lam;
Fren_x = abs( scalefactor * Gx.^2 / lambda_z );
Fren_y = abs( scalefactor * Gy.^2 / lambda_z );

if ( Fren_x <= Nx && Fren_y <= Ny ),
    % small Fresnel number case => far field
    %fprintf('Far Field\n');
    
    % set sample spacing of output beam
    if ~isempty(dxo) && ~isempty(dyo),
        if dxo > abs(lambda_z)/Gx || dyo > abs(lambda_z)/Gy,
            error('far field grid space too large');
        else
            dxout = dxo;
            dyout = dyo;
        end
    else % default values
        dxout = abs(lambda_z./Gx);
        dyout = abs(lambda_z./Gy);
    end

    % 
    beam = ApplyCurvature(beam,xin,yin,scalefactor/z,lam);
    beamout = Fraunhofer(beam,dx,dy,dxout,dyout,lambda_z);
        
    curvout = 1./z;
    xout = (-Nx/2:Nx/2-1)'*dxout;
    yout = (-Ny/2:Ny/2-1)'*dyout;
    
elseif ( Fren_x > Nx && Fren_y > Ny ),
    % large Fresnel number case => near field
    %fprintf('Near Field\n');
    
    % check scale of output grid to input grid, and apply geometric
    % curvature as necessary
    if ~isempty(dxo) && ~isempty(dyo),
        if abs( dxo/dx - dyo/dy ) <= 1e-6, % same scale
            scalefactor = sign(scalefactor) * dxo/dx;
            geomCurv    = ( scalefactor - 1 ) ./ z;
            beam = ApplyCurvature(beam, xin, yin, curv - geomCurv, lam);
            dxout = dxo;
            dyout = dyo;
        else
            error('Propagate: grid scale factor error');
        end
    else
        geomCurv = curv;
        dxout = abs(scalefactor) * dx;
        dyout = abs(scalefactor) * dy;
    end
    
    lambda_z_eff = lambda_z./scalefactor;
    
    % Fresnel propagation
    beamout = Fresnel;
    curvout = geomCurv./scalefactor;
    xout = (-Nx/2:Nx/2-1)'*dxout;
    yout = (-Ny/2:Ny/2-1)'*dyout;

else
    error('Propagate: Fresnel number conflict');
end
    
if applycurv,
    %fprintf('curvout = %f\n',curvout);
    beamout = ApplyCurvature(beamout, xout, yout, curvout, lam);
    curvout = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function beamout = Fraunhofer(beam,dx,dy,dxout,dyout,lambda_z)

        [Ny, Nx] = size(beam);
        Gx = Nx*dx; Gy = Ny*dy;

        ctmp = (-i./lambda_z) .* (dx*dy);
        dfx  = dxout./abs(lambda_z);
        dfy  = dyout./abs(lambda_z);

        if ( abs(Gx*dfx - 1) < 1e-9 ) && ( abs(Gy*dfy - 1) < 1e-9 ),
            beam = ifftshift(beam);
            beam = fft2(beam);
            beam = fftshift(beam);

        else
            beamin = beam;
            
            beam = czt2d(beam, dx, dy, dfx, dfy);

            beamintest = czt2d(beam, dfx, dfy, dx, dy, false);
            beamintestinverse = czt2d(beam, dfx, dfy, dx, dy, true);
            
            %save testfraunhofer.mat beamin beamintest beamintestinverse beam dx dy dfx dfy;
            
        end

        beamout = ctmp.*beam;

    end % Fraunhofer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function beamout = Fresnel %(beam,dx,dy,lambda_z_eff,scalefactor)
        
        %[Ny, Nx] = size(beam);
        %Gx = Nx*dx; Gy = Ny*dy;

        dfx = 1./Gx;
        dfy = 1./Gy;
        fx  = (-Nx/2:Nx/2-1)'*dfx;
        fy  = (-Ny/2:Ny/2-1)'*dfy;

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
        else
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
