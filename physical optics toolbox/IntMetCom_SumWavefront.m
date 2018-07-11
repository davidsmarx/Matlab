function status = IntMetCom_SumWavefront(h,ThreadID,BeamAID,BeamBID)
% status = IntMetCom_SumWavefront(h,ThreadID,BeamAID,BeamBID)
%
% cumulutive sum: BeamAID = BeamAID + BeamBID

% check if multi threading or not
if ThreadID == 0,
    status = h.SumWavefront(BeamAID,BeamBID);
else,
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('SumWavefront'),...
        [0 0],int32([0 0]),BeamAID,BeamBID,' ');
end

if status ~= 0,
    error(['SumWavefront error: ' num2str(status)]);
end
