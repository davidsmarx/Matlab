function status = IntMetCom_WavefrontPower(h,BeamID,ThreadID,powID)
% status = IntMetCom_WavefrontPower(h)
% returns current power in standard (single-threaded) wavefront
%
% status = IntMetCom_WavefrontPower(h,BeamID)
% returns current power in multi-threaded wavefront BeamID, assumes BeamID
% is not in a currently active thread.
%
% status = IntMetCom_WavefrontPower(h,BeamID,ThreadID,powID)
% adds a command to ThreadID to store the wavefront power of BeamID during
% thread execution in location powID
%

if ~exist('BeamID','var'), BeamID = 0; end
if ~exist('ThreadID','var'), ThreadID = 0; end

if ThreadID == 0,
    status = h.WavefrontPower(BeamID);
    
else,
    
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('WavefrontPower'),...
        [0 0],int32([powID 0]),BeamID,0,' ');
    if status ~= 0,
        error(['AddThreadCommand error: ' num2str(status)]);
    end

end
