function status = IntMetCom_CopyWavefront(h,ThreadID,BeamAID,BeamBID)
% status = IntMetCom_CopyWavefront(h,ThreadID,BeamAID,BeamBID)
%
% copies from BeamAID to BeamBID

% check if multi threading or not
if isempty(ThreadID) || ThreadID == 0,
    status = h.CopyWavefront(BeamAID,BeamBID);
else,
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('CopyWavefront'),...
        double([0 0]),int32([0 0]),BeamAID,BeamBID,' ');
end

if status ~= 0,
    error(['AddThreadCommand error: ' num2str(status)]);
end

% AddThreadCommand = int32 AddThreadCommand(handle, int32, int32,
% SafeArray(double), SafeArray(int32), int32, int32, string)
	