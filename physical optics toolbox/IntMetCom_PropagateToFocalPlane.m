function status = IntMetCom_PropagateToFocalPlane(...
    h,focuslens,dxout,dyout,ThreadID,BeamAID)
% status = IntMetCom_PropagateToFocalPlane(h,focuslens)
% status = IntMetCom_PropagateToFocalPlane(h,focuslens,dxout,dyout)
% status = IntMetCom_PropagateToFocalPlane(h,focuslens,dxout,dyout,ThreadID,BeamID)
% status = IntMetCom_PropagateToFocalPlane(h,focuslens,[],[],ThreadID,BeamID)
%
% focuslens is a struct:
%    focuslens.f = focal length
%    focuslens.D = clear aperture diameter

if exist('ThreadID','var') & ThreadID > 0,

    error('multi-thread version not yet implemented');
    
else,
    % single-threaded, default BeamID = 0
    if ~exist('BeamID','var'),
        BeamID = 0;
    end
    
    status = h.PropagateToFocalPlane(focuslens.f,focuslens.D,dxout,dyout,...
        true,BeamID);

end

if status ~= 0,
    error(['PropagateToFocalPlane error: ' num2str(status)]);
end
