function IntMetCom_PutField(varargin)
% IntMetCom_PutField(h,field,dx,dy,wavelength,curv)
% IntMetCom_PutField(h,field,dx,dy,wavelength,curv,BeamID)
% IntMetCom_PutField(h, S, BeamID)
%
% dx, dy, wavelength, curv are optional, default values are:
%    dx = 15*MM/nx
%    dy = 15*MM/ny
%    wavelength = 1.319*UM
%    curv = 0
% where [ny,nx] = size(field)
%
% there is no multi-threaded version of PutField

unitsdefinitions;

switch nargin,
    % h,field,dx,dy,wavelength,curv,BeamID)
    case 2,
        [h, S] = deal(varargin{:});
        [field, dx, dy, wavelength, curv] =...
            deal(S.field, S.dx, S.dy, S.wavelength, S.curv);
        BeamID = 0;
        
    case 3,
        [h, S, BeamID] = deal(varargin{:});
        [field, dx, dy, wavelengh, curv] =...
            deal(S.field, S.dx, S.dy, S.wavelength, S.curv);
        
    case 6,
        [h, field, dx, dy, wavelength, curv] = deal(varargin{:});
        BeamID = 0;
        
    case 7,
        [h, field, dx, dy, wavelength, curv, BeamID] = deal(varargin{:});
        
    otherwise,
        error('usage');
end

[ny, nx] = size(field);

if ~exist('dx','var') || isempty(dx), dx = 15*MM./nx; end
if ~exist('dy','var') || isempty(dy), dy = 15*MM./ny; end
if ~exist('wavelength','var') || isempty(wavelength), wavelength = 1.319*UM; end
if ~exist('curv','var') || isempty(curv), curv = 0; end

status = h.wavefront_put(real(field),imag(field),...
    dx,dy,wavelength,curv,BeamID);

IntMetComCheckStatus(status);

