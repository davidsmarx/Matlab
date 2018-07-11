function status = IntMetCom_WriteWavefrontUNF(h,filename,ThreadID,BeamID)
% status = IntMetCom_WriteWavefrontUNF(h,filename)
% status = IntMetCom_WriteWavefrontUNF(h,filename,ThreadID,BeamID)

if ~exist('ThreadID','var') | isempty(ThreadID),
    ThreadID = 0;
end

if ~exist('BeamID','var'),
    BeamID = 0;
end

if ThreadID == 0,
    status = h.WriteWavefrontUNF(filename,BeamID);

else,
    status = h.AddThreadCommand(ThreadID,IntMetComCMD('WriteWavefrontUNF'),...
        [0 0],int32([0 0]),BeamID,0,filename);    
end

IntMetComCheckStatus(status);
