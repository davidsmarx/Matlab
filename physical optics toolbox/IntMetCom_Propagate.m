function status = IntMetCom_Propagate(h,...
    distance,ThreadID,BeamID,dxout,dyout,applycurv)
% status = IntMetCom_Propagate(h,distance)
% status = IntMetCom_Propagate(h,distance,ThreadID,BeamID)
% status =...
%     IntMetCom_Propagate(h,distance,[],[],dxout,dyout,applycurv)
% status =...
%     IntMetCom_Propagate(h,distance,ThreadID,BeamID,dxout,dyout,applycurv)
%
% default dxout, dyout = wavefront's current dx, dy
% default applycurv = true

if exist('ThreadID','var') & ThreadID > 0,
    % multi-threaded
    if ~exist('dxout','var') | ~exist('dyout','var')
        rp = [distance 0];
    else
        rp = [distance dxout dyout];
    end
    if ~exist('applycurv','var'),
        ip = int32([1 0]); % default = true
    else
        ip = int32([applycurv 0]);
    end
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('Propagate'),...
        rp, ip, BeamID, 0, ' ');

else,
    % single-threaded, default BeamID = 0
    if ~exist('BeamID','var') || isempty(BeamID),
        BeamID = 0;
    end
    if ~exist('dxout','var') || ~exist('dyout','var')
        [nx, ny, dxout, dyout, curv, lam] = h.GetWavefrontParms(BeamID);
    end
    if ~exist('applycurv','var'),
        applycurv = true; % default = true
    end

    status = h.PropagateExt(distance,dxout,dyout,applycurv,BeamID);

end

if status ~= 0,
    error(['Propagate error: ' num2str(status)]);
end
