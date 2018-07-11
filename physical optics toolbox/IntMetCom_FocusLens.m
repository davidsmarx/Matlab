function status = IntMetCom_FocusLens(h,focuslens,dxout,dyout,ThreadID,BeamAID)
% status = IntMetCom_FocusLens(h,focuslens,dxout,dyout)
% status = IntMetCom_FocusLens(h,focuslens,dxout,dyout,ThreadID,BeamAID)
%
% focuslens is a struct:
%    focuslens.f = focal length
%    focuslens.D = clear aperture diameter

% !!call FocusLens( ObjectData%wvfront, focuslens_f, focuslens_D, dxout, dyout )
% call FocusLens( cmdCurrent%pwA, cmdCurrent%rP(1,1), cmdCurrent%rP(1,2),&
% cmdCurrent%rP(1,3), cmdCurrent%rP(1,4) )

if exist('ThreadID','var') & ~isempty(ThreadID) & ThreadID > 0,

    rP = [focuslens.f focuslens.D dxout dyout];
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('FocusLens'),...
        rP,int32([0 0]),BeamAID,0,' ');

else, 
    % single-threaded
    if ~exist('BeamAID','var'),
        BeamAID = 0;
    end
    status = h.focuslens(focuslens.f,focuslens.D,dxout,dyout,BeamAID);
end

if status ~= 0,
    error(['FocusLens error: ' num2str(status)]);
end
