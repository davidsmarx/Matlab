function phi = phaseunwrap2d(psi, varargin)
% phi = phaseunwrap2d(psi, sOptions)
%
% unweighted MSE direct solution to the 2d phase unwrapping problem
% maybe add some other types of solutions in the future
%
% input psi = wrapped phase map array
%
% output phi = unwrapped phase
%
% (optional) sOptions = struct(...
%                   'RetainOriginalPhase', true... (default) replace
%                               output phase with input psi +/ integer 2*pi
%                   );

sOptions = ValidateOptions(varargin{:});

[Ny, Nx] = size(psi);

% direct solution to equation 5.31 in Ghiglia & Pritt
% expand to N+1 x N+1 before differencing back to NxN

% the diagonals, (:) is row,column order
aij = -4*ones(Nx*Ny,1);
aipj = [ones(Ny,Nx-1) zeros(Ny,1)]; 
aimj = [zeros(Ny,1) ones(Ny,Nx-1)];
aijp = [ones(Ny-1,Nx); zeros(1,Nx)];
aijm = [zeros(1,Nx); ones(Ny-1,Nx)];

% build sparse matrix
B = -[circshift(aimj(:),-Ny) circshift(aijm(:),-1) aij(:) circshift(aijp(:),1) circshift(aipj(:),Ny)];
A = spdiags(B, [-Ny -1 0 1 Ny], Nx*Ny, Nx*Ny);

% % Neumann boundary conditions = 1st derivative = 0 across boundary
% % index into sparse array for ii,jj with ii,jj origin at center of image
p = @(ii,jj) ((ii+Nx/2)*Nx + jj+Ny/2+1);
% down the left boundary: A(1:N,:), d/dx => phi(ii+1,jj)-phi(ii,jj)
for jj = -Ny/2 : Ny/2-1,
    ii = -Nx/2;
    
    A(p(ii,jj),p(ii,jj)) = A(p(ii,jj),p(ii,jj)) + 1;
end

% up the right boundary: A(end-N+1:end,:)
for jj = -Ny/2 : Ny/2-1  % p = N.^2-N+1:N.^2
    ii = Nx/2-1;
    
    A(p(ii,jj),p(ii,jj)) = A(p(ii,jj),p(ii,jj)) + 1;
end

% across top boundary: A(1:N:end,:)
for ii = -Nx/2 : Nx/2-1,
    jj = -Ny/2;

    A(p(ii,jj),p(ii,jj)) = A(p(ii,jj),p(ii,jj)) + 1;
end

% left-to-right across bottom boundary: A(N:N:end,:)
for ii = -Nx/2 : Nx/2-1,
    jj = Ny/2-1;

    A(p(ii,jj),p(ii,jj)) = A(p(ii,jj),p(ii,jj)) + 1;
end
%%%%%%%%%%%%%%%%%%%%% end Neumann boundary conditions

% rhs
del_x_ij  = mod2pi(diff([psi psi(:,end)],1,2));
del_x_imj = mod2pi(diff([psi(:,1) psi],1,2));
del_y_ij  = mod2pi(diff([psi; psi(end,:)]));
del_y_ijm = mod2pi(diff([psi(1,:); psi]));

rhs = (del_x_ij - del_x_imj) + (del_y_ij  - del_y_ijm);
rhs = -rhs(:);

% 
% icoluse = ~isnan(rhs);
% Ap = A(:,icoluse);
% size(Ap)
rhs(isnan(rhs)) = 0;
% phip = Ap \ rhs;
% size(phip)
% phipp = zeros(size(rhs));
% phipp(icoluse) = phip;
% phipp = reshape(phipp,Ny,Nx);
% figure,imagesc(phipp), axis image

% solve
phi = A \ rhs;

phi = reshape(phi,Ny,Nx);

if isfield(sOptions, 'RetainOriginalPhase') && sOptions.RetainOriginalPhase,
    % retain original phase values, just +/- 2pi where indicated
    pixplus = phi - psi > pi;
    pixminu = phi - psi < -pi;
    phi_unw = psi;
    phi_unw(pixplus) = psi(pixplus) + 2*pi;
    phi_unw(pixminu) = psi(pixminu) - 2*pi;
   
    phi = phi_unw;
end

end % main

function sOptions = ValidateOptions(varargin)

    % defaults
    sOptions = struct(...
        'RetainOriginalPhase', true ...
        );
    
    if nargin > 0 && isstruct(varargin{1}),
        sIn = varargin{1};
        fnames = fieldnames(sIn);
        for ifname = 1:length(fnames),
            sOptions.(fnames{ifname}) = sIn.(fnames{ifname});
        end
        
    end
    
end % ValidateOptions