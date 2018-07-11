function status = IntMetCom_ThinLens(h,fl,ThreadID,BeamID)
% status = IntMetCom_ThinLens(h,focuslens)
% status = IntMetCom_ThinLens(h,focuslens,ThreadID,BeamID)
%
% focuslens is a struct:
%    focuslens.f = focal length
%    focuslens.D = clear aperture diameter
%    focuslens.xc, focuslens.yc (optional) = lens decenter

if ~isfield(fl,'xc'), fl.xc = 0; end
if ~isfield(fl,'yc'), fl.yc = 0; end

if exist('ThreadID','var') & ThreadID > 0,
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('ThinLens'),...
        [fl.f, fl.D, fl.xc, fl.yc],int32([0 0]),BeamID,0,' ');

else,
    % single-threaded
    if ~exist('BeamID','var'),
        BeamID = 0;
    end
    
    status = h.ThinLens(fl.f,fl.D,fl.xc,fl.yc,BeamID);
end

if status ~= 0,
    error(['Propagate error: ' num2str(status)]);
end
