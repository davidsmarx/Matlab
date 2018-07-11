function status = IntMetCom_CopyWavefrontAmp(h,ThreadID,BeamAID,BeamBID)
% status = IntMetCom_CopyWavefrontAmp(h,ThreadID,BeamAmpID,BeamPhaseID)
%
% copies the amplitude from BeamAmpID to BeamPhaseID, so that BeamPhaseID
% is overwritten with the phase of BeamPhaseID and the amplitude of
% BeamAmpID. BeamAmpID is not overwritten

% check if multi threading or not
if isempty(ThreadID) | ThreadID == 0,
    status = h.CopyWavefrontAmp(BeamAID,BeamBID);
else,
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('CopyWavefrontAmp'),...
        [0 0],int32([0 0]),BeamAID,BeamBID,' ');
end

IntMetComCheckStatus(status);