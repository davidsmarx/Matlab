function [ A, dx, dy, lambda, curv, x, y ] = ReadWavefrontUNF( fileName, fun_handle, varargin )
%===============================================================================
%
%	[ A, dx, dy, lambda, curv, x, y ] = ReadWavefrontUNF( fileName )
%
%	Reads an unformatted file containing wavefront data...probably from a
%	Fortran diffraction simulation (fileName would be of the form "name.unf").
%	The file would have been written by WriteWavefrontUNF.F95.
%
%	A      = complex-valued, Ny by Nx wavefront amplitude array
%	dx, dy = Sample spacing
%	lambda = Wavelength
%	curv   = 1/wavefront curvature
%	x, y   = Coordinate values of sample points (column vectors)
%
%				Units: A      - sqrt(watts)/meter ( A.*conjg(A) = intensity )
%						 dx, dy - meters
%						 lambda - meters
%						 curv   - 1/meters
%						 x, y   - meters
%
%		22 June '01
%
%   ReadWavefrontUNF(fileName)
%   ReadWavefrontUNF(fileName,fun_handle)
%   ReadWavefrontUNF(fileName,fun_handle,options)
%   ReadWavefrontUNF(fileName,[],options)
%
%   with no output arguments, make a figure showing the wavefront. options
%   are options to imagesc(). What is plotted is fun_handle(wavefront). For
%   example, fun_handle = @abs (default).
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

fileName = CheckInput(fileName);

[ fid, message ] = fopen( fileName, 'r', 'native' );
if ( fid == -1 )
   error( ['failed to open file ' fileName] )
end

try,
    Nsamp = fread( fid, 2, 'int32' );			% Number of samples in beam grid
    Nx = Nsamp(1);
    Ny = Nsamp(2);

    lambda = fread( fid, 1, 'double' );			% Wavelength
    curv   = fread( fid, 1, 'double' );			% Wavefront curvature

    x = fread( fid, Nx, 'double' );				% x-coordinate values
    y = fread( fid, Ny, 'double' );				% y-coordinates

    dx = ( x(Nx) - x(1) )/(Nx-1);
    dy = ( y(Ny) - y(1) )/(Ny-1);

    A = fread( fid, [ Ny, Nx ], 'double' );	% Real part of wave amplitude
    B = fread( fid, [ Ny, Nx ], 'double' );

catch, % make sure to close file if error
    status = fclose(fid);
    rethrow(lasterror);
end % try-catch

% A = A + i*fread( fid, [ Ny, Nx ], 'double' );
% memory saving version:
for ii = 1:Nx,
    A(:,ii) = A(:,ii) + i*B(:,ii);
end

status = fclose( fid );
if ( status ~= 0 )
   error( 'ReadWavefrontUNF: error closing file' )
end

% if no output arguments, show the wavefront in a figure
if nargout == 0,
    if ~exist('fun_handle','var') | isempty(fun_handle),
        fun_handle = @abs;
    end
    ShowWavefront(x,y,fun_handle(A),fileName,varargin{:});
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileName = CheckInput(fn)

if ~strcmp(fn(end-3:end),'.unf'),
    fileName = [fn '.unf'];
else,
    fileName = fn;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowWavefront(x,y,A,titlestr,varargin)

global MM;

imagesc(x/MM,y/MM,A,varargin{:}), axis image, colorbar,
title(pwd2titlestr(titlestr));
