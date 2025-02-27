function status = IntMetCom_CreateGausSource(h,beamwaistdiam,varargin)
% status = IntMetCom_CreateGausSource(h,beamwaistdiam)
% status = IntMetCom_CreateGausSource(h,beamwaistdiam,s_parms)
% status = IntMetCom_CreateGausSource(h,beamwaistdiam,s_parms,ThreadID,beamID)
% status = IntMetCom_CreateGausSource(h,beamwaistdiam,[],ThreadID,BeamID)
%
% h = handle to com server
% beamwaistdiam = beam waist diameter
% s_parms = optional struct with the following optional fields and
% defaults:
%     nx = x-dimension of wavefront matrix (default = 512)
%     ny = y-dimension (default = 512)
%     Gx = grid size in x (default = 15*MM), pixel spacing = Gx/nx
%     Gy = grid size in y (default = 15*MM), pixel spacing = Gy/ny
%     wavelength = wavelength (default = 1.319*UM)
%
% Wavefront Power = 1.0, by default in the COM method CreateGausSource
%
% when ThreadID and beamID are supplied, AddThreadCommand is used to add
% the CreateGausSource command to the thread.

[nx, ny, dx, dy, wavelength, ThreadID, BeamID] = ValidateInput(varargin{:});

if ThreadID == 0,
    status = h.CreateGausSource(beamwaistdiam,nx,ny,dx,dy,wavelength,BeamID);
else,
	status = h.AddThreadCommand(ThreadID,IntMetComCMD('CreateGausSource'),...
        [dx dy wavelength beamwaistdiam], int32([nx ny]),BeamID,0,'a');
end

IntMetComCheckStatus(status);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nx, ny, dx, dy, wavelength, ThreadID, BeamID] = ValidateInput(varargin)
% varargin{:} can be:
%   {s_parms}
%   {s_parms, ThreadID, beamID}
%   {[], ThreadID, beamID}

global MM UM;

% define defaults:
nx = 512; ny = 512;
dx = 15*MM/nx; dy = dx;
wavelength = 1.319*UM;
ThreadID = 0;
BeamID = 0;

if length(varargin) > 0,
    stmp = varargin{1};
    if isstruct(stmp),
        if isfield(stmp,'nx'), nx = stmp.nx; end
        if isfield(stmp,'ny'), ny = stmp.ny; end
        if isfield(stmp,'Gx'), dx = stmp.Gx/nx; end
        if isfield(stmp,'Gy'), dy = stmp.Gy/ny; end
        if isfield(stmp,'wavelength'), wavelength = stmp.wavelength; end
    end
end
if length(varargin) == 3,
    if ~isempty( varargin{2} ), ThreadID = varargin{2}; end
    if ~isempty( varargin{3} ), BeamID = varargin{3}; end
    if ThreadID < 0,
        error(['Invalid ThreadID: ' num2str(ThreadID)]);
    end
    if BeamID < 0,
        error(['Invalid BeamID: ' num2str(BeamID)]);
    end
end
