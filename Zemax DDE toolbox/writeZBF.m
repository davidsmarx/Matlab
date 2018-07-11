function emsg = writeZBF(filename, Ex, Ey, S)
% emsg = writeZBF(filename, Ex, Ey, S)
% Ex = complex field
% Ey = complex field (=0 or empty if unpolarized)
% S is a struct with fields:
%   Nx, Ny = ize of Ex and Ey
%   polflag = 0 => unpolarized, Ey = 0, 1 => polarized
%   dx, dy = ample pacing
%   zz = poition relative to pilot beam wait
%   zr = Rayleigh ditance for the pilot beam
%   ww = beam wait of the pilot beam
%   lam = wavelength [len unit]

U = CConstants;

emsg = '';

fid = fopen(filename, 'wb');
if fid == -1, error(['error opening ' filename]); end

% write first 9 parameter integer
iUnits = 0; UU = U.MM; % units
A = [1 S.nx S.ny S.polflag iUnits 0 0 0 0];
cnt = fwrite(fid, A, 'int32');
if cnt ~= 9, emsg = 'error writing integer parameters'; return, end

% write next 20 parameters double
A = [S.dx/UU S.dy/UU S.zwx/UU S.zrx/UU S.wwx/UU S.zwy/UU S.zry/UU S.wwy/UU S.wavelength/UU S.ior zeros(1, 10)];
cnt = fwrite(fid, A, 'double');
if cnt ~= 20, emsg = 'error writing double parameters'; return, end

% write the real and imag parts of the field
% need to transpose the field arrays to conform to zemax order
Ex = Ex.';
A = [real(Ex(:))'; imag(Ex(:))'];
cnt = fwrite(fid, A(:), 'double');

if S.polflag,
    Ey = Ey.';
    A = [real(Ey(:))'; imag(Ey(:))'];
    cnt = fwrite(fid, A(:), 'double');
end

fclose(fid);

return