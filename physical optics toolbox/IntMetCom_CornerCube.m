function structCcube = IntMetCom_CornerCube(h,Ccube,ThreadID,WavefrontID,direction,sidtilt)
% structCcube = IntMetCom_CornerCube(h,Ccube)
% structCcube = IntMetCom_CornerCube(h,Ccube,ThreadID,WavefrontID)
% structCcube = IntMetCom_CornerCube(h,Ccube,ThreadID,WavefrontID,direction,sidtilt)
% structCcube = IntMetCom_CornerCube(h,Ccube,[],[],direction,sidtilt)
% 
% structCcube = IntMetCom_CornerCube;
%    returns a structure with all the fields set to default values
%
% IntMetCom_CornerCube(h,Ccube) invokes the CornerCube method in the
% interface h, using the corner cube parameters stored in the struct Ccube
%
% structCcube = IntMetCom_CornerCube(h,Ccube,direction,sidtilt)
% direction = 'ns' or 'ew'
% sidtilt = struct with sidtilt.x and sidtilt.y = rotations about x and y
%    axes
% corner_cube_orientations is called to get the rotation matrix for CCube
%
% The fields in the Ccube structure are:
%    size
%    shape ('square' (default), 'triangle', 'sim')
%    edgelen(1:3) = relative length of each roofline. (default = [1 1 1]).
%         The length of each roofline will be size*edgelen(ii).
%    xc  = vertex offset in x (default = 0)
%    yc  = vertex offset in y (default = 0)
%    spin = clocking of corner cube (default = 0)
%    gapwidth(1:3) = size of gap along each roofline (default = 0)
%    dihedral(1:3) = dihedral (rad) error of each roofline
%    rotmatrix(1:3,1:3) = rotation matrix input to Diffraction code.
%       default = eye(3);
%    edgedir(1:3,1:3) = matrix of three column vectors, where each vector is a
%       direction vector indicating the direction of a roof line.
%
% ThreadID = Corner Cube object is added to thread ThreadID. ThreadID = 0
% for single threaded (default)
% WavefrontID = Corner Cube object applied to this wavefront.
%
% ThreadID and WavefrontID are positive integers.

unitsdefinitions;

% if no input arguments, return default Ccube struct
if nargin == 0,
    % return default struct if no input args
    structCcube = struct(...
        'size',75*MM/7,...
        'edgelen',[1 1 1],...
        'shape','square',...
        'xc',0,...
        'yc',0,...
        'spin',0,...
        'gapwidth',[0 0 0]',...
        'dihedral',[0 0 0]',...
        'rotmatrix',eye(3),...
        'edgedir',eye(3) );
    return;
end

% set defaults
if ~exist('ThreadID','var') | isempty(ThreadID),
    ThreadID = 0; % single threaded
end
if ~exist('WavefrontID','var') | isempty(WavefrontID),
    WavefrontID = 0;
end

% if direction and sidtilt are input, orient the corner cube
if exist('direction','var')
    switch lower(direction)
        case {'ns','c4'}
            Ccube.rotmatrix = corner_cube_orientations(4,sidtilt.x,sidtilt.y);
        case {'ew','c3'}
            Ccube.rotmatrix = corner_cube_orientations(3,sidtilt.x,sidtilt.y);
        otherwise,
            error(['unknown direction: ' direction]);
    end
end

% if ThreadID == 0, run single threaded
if ThreadID == 0,
    status = InvokeCornerCube(h,Ccube,WavefrontID);
else,
    status = InvokeCornerCubeMT(h,Ccube,ThreadID,WavefrontID);
end

IntMetComCheckStatus(status);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = InvokeCornerCube(h,Ccube,WavefrontID)

status = h.CornerCube(Ccube.size,Ccube.shape,...
    Ccube.xc,Ccube.yc,Ccube.spin,...
    Ccube.gapwidth,Ccube.rotmatrix,Ccube.dihedral,Ccube.edgelen,...
    WavefrontID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = InvokeCornerCubeMT(h,Ccube,ThreadID,WavefrontID)

status = h.AddThreadCommand(ThreadID,IntMetComCMD('CornerCube'),...
    [
    Ccube.size Ccube.xc Ccube.yc
    Ccube.spin 0 0
    Ccube.gapwidth(:)'
    Ccube.dihedral(:)'
    Ccube.edgelen(:)'
    Ccube.rotmatrix
    ], int32([0 0]), WavefrontID, 0, Ccube.shape);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function flagMT = ValidateID(ThreadID,WavefrontID)
% 
% flagMT = 0;
% if ~isempty(ThreadID) & ~isempty(WavefrontID)
%     if ThreadID > 0 & WavefrontID > 0,
%         flagMT = 1;
%     end
% end
