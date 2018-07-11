function [ Ar, dr, lambda, r, intensity] = ReadWavefrontUNFradial( fileName )
% [ Ar, dr, lambda, r, intensity] = ReadWavefrontUNFradial( fileName )

[ A, dx, dy, lambda, curv, x, y ] = ReadWavefrontUNF( fileName );

if curv ~= 0, error('need to apply curvature'); end

r = y(y>=0);
Ar = A(y>=0,x==0);
dr = dy;

intensity = wavefront_intensity(A,dx,dy);