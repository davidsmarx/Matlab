function [Aout, phi] = phase_retrieval_TIE(I1, I2, dx, dy, wavelength, zprop, sOptions)
% [Aout, Phiout] = phase_retrieval_TIE(I1, I2, dx, dy, wavelength, zprop, options)
%
% I1 and I2 are intensity arrays = amp^2
% the wavefront propagates from I1 to I2
% zprop = distance from I1 to I2
%
% references:
%     p. 181 of notebook (1/12/09)
%     T.E. Gureyev, "Composite Techniques for Phase Retrieval in the
%        Fresnel Region," Optics Comm., vol. 220, p. 49, 2003.
%     Numerical Recipes in C chapter on PDE's
%
% options is a struct
%     'debug' (default = off)

NEUMANN = false;

if dx ~= dy,
    error('dx ~= dy not implemented');
end
[Ny, Nx] = size(I1);
if any(size(I2) ~= [Ny, Nx]),
    error('input I1 and I2 must be same size');
end
if Ny ~= Nx
    error('input I1 and i2 must be square');
end
N = Nx;

if ~exist('sOptions','var')
    sOptions = struct;
end
sOptions = CheckOptions(sOptions);

%TIE
% ii,jj is column,row, but matlab (:) operator is row,column
I0ipjtmp = 0.5*(I1(:,2:end) + I1(:,1:end-1));
I0ipj = [I0ipjtmp I0ipjtmp(:,end)];
I0imj = [I0ipjtmp(:,1) I0ipjtmp];
I0ijptmp = 0.5*(I1(2:end,:) + I1(1:end-1,:));
I0ijp = [I0ijptmp; I0ijptmp(end,:)];
I0ijm = [I0ijptmp(1,:); I0ijptmp];

aipj = [I0ipjtmp zeros(N,1)];
aimj = [zeros(N,1) I0ipjtmp];
aijp = [I0ijptmp; zeros(1,N)];
aijm = [zeros(1,N); I0ijptmp];

% the diagonals, (:) is row,column order
aij = -(I0ipj(:) + I0imj(:) + I0ijp(:) + I0ijm(:));
aipj = aipj(:);
aimj = aimj(:);
aijp = aijp(:);
aijm = aijm(:);

% build sparse matrix
B = -[circshift(aimj,-N) circshift(aijm,-1) aij circshift(aijp,1) circshift(aipj,N)];
%B = -[aimj aijm aij aijp aipj];
S = spdiags(B, [-N -1 0 1 N], N^2, N^2);

% rhs of pde equation
rhs = (2*pi*dx*dx/(wavelength*zprop))*(I2(:) - I1(:));

% index into sparse array for ii,jj with ii,jj origin at center of image
p = @(ii,jj) ((ii+Nx/2)*Nx + jj+Ny/2+1);

% Boundary Conditions
if NEUMANN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Neumann boundary conditions
    % % boundary conditions on the outward gradient of phi
    % % see p. 189 of notebook and Gureyev 1995.

    % down the left boundary: S(1:N,:), d/dx => phi(ii+1,jj)-phi(ii,jj)
    for jj = -Ny/2+1 : Ny/2-2,
        ii = -Nx/2;

        sin_over_cos = jj/ii;

        S(p(ii,jj),p(ii+1,jj)) = 2*S(p(ii,jj),p(ii+1,jj));
        %     S(p(ii,jj),p(ii,jj+1)) = (1+sin_over_cos)*S(p(ii,jj),p(ii,jj+1));
        %     S(p(ii,jj),p(ii,jj-1)) = (1-sin_over_cos)*S(p(ii,jj),p(ii,jj-1));
    end

    % up the right boundary: S(end-N+1:end,:)
    for jj = -Ny/2+1 : Ny/2-2  % p = N.^2-N+1:N.^2
        ii = Nx/2-1;

        sin_over_cos = jj/ii;

        S(p(ii,jj),p(ii-1,jj)) = 2*S(p(ii,jj),p(ii-1,jj));
        %     S(p(ii,jj),p(ii,jj+1)) = (1-sin_over_cos)*S(p(ii,jj),p(ii,jj+1));
        %     S(p(ii,jj),p(ii,jj-1)) = (1+sin_over_cos)*S(p(ii,jj),p(ii,jj-1));
    end

    % across top boundary: S(1:N:end,:)
    for ii = -Nx/2+1 : Nx/2-2,
        jj = -Ny/2;

        cos_over_sin = ii/jj;

        %     S(p(ii,jj),p(ii-1,jj)) = (1-cos_over_sin)*S(p(ii,jj),p(ii-1,jj));
        %     S(p(ii,jj),p(ii+1,jj)) = (1+cos_over_sin)*S(p(ii,jj),p(ii+1,jj));

        S(p(ii,jj),p(ii,jj+1)) = 2*S(p(ii,jj),p(ii,jj+1));

    end

    % left-to-right across bottom boundary: S(N:N:end,:)
    for ii = -Nx/2+1 : Nx/2-2,
        jj = Ny/2-1;

        cos_over_sin = ii/jj;

        %     S(p(ii,jj),p(ii-1,jj)) = (1+cos_over_sin)*S(p(ii,jj),p(ii-1,jj));
        %     S(p(ii,jj),p(ii+1,jj)) = (1-cos_over_sin)*S(p(ii,jj),p(ii+1,jj));

        S(p(ii,jj),p(ii,jj-1)) = 2*S(p(ii,jj),p(ii,jj-1));

    end

    % the corners (assumes square grid => Nx = Ny, otherwise need to account
    % for sin ~= cos at the corners
    ii = -Nx/2; jj = -Nx/2;
    S(p(ii,jj),p(ii+1,jj)) = 2*S(p(ii,jj),p(ii+1,jj));
    S(p(ii,jj),p(ii,jj+1)) = 2*S(p(ii,jj),p(ii,jj+1));
    ii =  Nx/2-1; jj = -Nx/2;
    S(p(ii,jj),p(ii-1,jj)) = 2*S(p(ii,jj),p(ii-1,jj));
    S(p(ii,jj),p(ii,jj+1)) = 2*S(p(ii,jj),p(ii,jj+1));
    ii = -Nx/2; jj = Nx/2-1;
    S(p(ii,jj),p(ii+1,jj)) = 2*S(p(ii,jj),p(ii+1,jj));
    S(p(ii,jj),p(ii,jj-1)) = 2*S(p(ii,jj),p(ii,jj-1));
    ii = Nx/2-1; jj = Nx/2-1;
    S(p(ii,jj),p(ii-1,jj)) = 2*S(p(ii,jj),p(ii-1,jj));
    S(p(ii,jj),p(ii,jj-1)) = 2*S(p(ii,jj),p(ii,jj-1));

    %%%%%%%%%%%%%%%%%%%%% end Neumann boundary conditions

else % dirichlet
    %%%%%%%%%%%%%%%%% Dirichlet boundary conditions
    for jj = -Ny/2 : Ny/2-1,
        ii = -Nx/2;

        S(p(ii,jj),:) = 0;
        S(p(ii,jj),p(ii,jj)) = 1;
        rhs(p(ii,jj)) = 0;
    end

    % up the right boundary: S(end-N+1:end,:)
    for jj = -Ny/2 : Ny/2-1  % p = N.^2-N+1:N.^2
        ii = Nx/2-1;

        S(p(ii,jj),:) = 0;
        S(p(ii,jj),p(ii,jj)) = 1;
        rhs(p(ii,jj)) = 0;
    end

    % right-to-left across top boundary: S(1:N:end,:)
    for ii = -Nx/2 : Nx/2-1,
        jj = -Ny/2;

        S(p(ii,jj),:) = 0;
        S(p(ii,jj),p(ii,jj)) = 1;
        rhs(p(ii,jj)) = 0;
    end

    % left-to-right across bottom boundary: S(N:N:end,:)
    for ii = -Nx/2 : Nx/2-1,
        jj = Ny/2-1;

        S(p(ii,jj),:) = 0;
        S(p(ii,jj),p(ii,jj)) = 1;
        rhs(p(ii,jj)) = 0;

    end
    %%%%%%%%%%%%%%%%%% end Dirichlet boundary conditions

end % boundary conditions

% need a constant
% S(p(-Nx/2,  -Ny/2),:) = 0;   S(p(-Nx/2,  -Ny/2),  p(-Nx/2,  -Ny/2)) = 1;   rhs(p(-Nx/2,  -Ny/2)) = 0;
% S(p( Nx/2-1,-Ny/2),:) = 0;   S(p( Nx/2-1,-Ny/2),  p( Nx/2-1,-Ny/2)) = 1;   rhs(p( Nx/2-1,-Ny/2)) = 0;
% S(p(-Nx/2,   Ny/2-1),:) = 0; S(p(-Nx/2,   Ny/2-1),p(-Nx/2,   Ny/2-1)) = 1; rhs(p(-Nx/2,   Ny/2-1)) = 0;
% S(p( Nx/2-1, Ny/2-1),:) = 0; S(p( Nx/2-1, Ny/2-1),p( Nx/2-1, Ny/2-1)) = 1; rhs(p( Nx/2-1, Ny/2-1)) = 0;
% S(1,1) = 0;
% rhs(1) = 0;

% solve
%phi = ( S \ rhs );
phi = bicgstab(S,rhs,1e-4,100);

% offset by arbitrary constant
% also remove tilt
x = (-Nx/2:Nx/2-1);
y = (-Ny/2:Ny/2-1);
[X, Y] = meshgrid(x,y);
[Z, phi_residual, rmsresidual, R] = zernikefit(X(:), Y(:), phi, [1 2 3]);

% return to image array format
phi = reshape(phi_residual,Ny,Nx);

% create a complex field with input amplitude and calculated phase
Aout = sqrt(I1).*exp(j*phi);

%
% if strcmp(sOptions.debug,'on')
%     x = (-Nx/2:Nx/2-1)*dx;
%     y = (-Ny/2:Ny/2-1)*dy;
%     MM = 1e-3;
%     figure('position',[847 19 560 420]), imagesc(x/MM, y/MM, angle(Aout)/pi),
%     axis image,
%     colorbar, colorbartitle('Phase [\pi rad]')
%     title('TIE Calculated Phase')
% end

end % main

function sOptions = CheckOptions(sIn)
% default values:
sOptions = struct(...
    'debug','off'...
    );

if ~isstruct(sIn),
    return
end

if isfield(sIn,'debug'),
    sOptions.debug = sIn.debug;
end
% add more options

end % CheckOptions
