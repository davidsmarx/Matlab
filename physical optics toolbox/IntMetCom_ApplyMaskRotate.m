function status = IntMetCom_ApplyMaskRotate(h,bsmask,masktheta,xc,yc,direction,...
    ThreadID,BeamID)
% status = IntMetCom_ApplyMaskRotate(h,bsmask,masktheta,xc,yc,direction)
% status = IntMetCom_ApplyMaskRotate(h,bsmask,masktheta,xc,yc,direction,...
%    ThreadID,BeamID)

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
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('ApplyMaskRotate'),...
        [bsmask(:)' masktheta xc yc],int32([0 0]),BeamID,0,dd);

else,
    % single-threaded
    if ~exist('BeamID','var'),
        BeamID = 0;
    end
    status = h.ApplyMaskRotate(bsmask(1),bsmask(2),bsmask(3),masktheta,...
        xc,yc,dd,BeamID);
    
end

IntMetComCheckStatus(status);
