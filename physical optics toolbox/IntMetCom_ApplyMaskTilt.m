function status = IntMetCom_ApplyMaskTilt(h,bsmask,masktheta,xc,yc,tilt,direction,...
    ThreadID,BeamID)
% status = IntMetCom_ApplyMaskTilt(h,bsmask,masktheta,xc,yc,tilt,direction)
% status = IntMetCom_ApplyMaskTilt(h,bsmask,masktheta,xc,yc,tilt,direction,...
%    ThreadID,BeamID)
%
% bsmask = [length, width, offset]
% masktheta = rotation about propagation axis
% xc, yc = translation offset of mask
% tilt = [angle, orientation] of tilt of mask with respect to propagation
%    axis. orientation is angle between tilt axis and x-axis

% translate various direction options
switch lower(direction)
    case {'ns','vv',2},
        dd = 'ns';
    case {'ew','hh',1},
        dd = 'ew';
    otherwise,
        error(['unknown direction: ' direction]);
end

% if a valid ThreadID, use multi-threading call
if exist('ThreadID','var') & ~isempty(ThreadID) & ThreadID > 0,
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('ApplyMaskTilt'),...
        [bsmask(:)' masktheta xc yc tilt(:)'],int32([0 0]),BeamID,0,dd);

else,
    % single-threaded
    if ~exist('BeamID','var'),
        BeamID = 0;
    end
    status = h.ApplyMaskTilt(bsmask(1),bsmask(2),bsmask(3),masktheta,...
        xc,yc,tilt(1),tilt(2),dd,BeamID);
    
end

IntMetComCheckStatus(status);
