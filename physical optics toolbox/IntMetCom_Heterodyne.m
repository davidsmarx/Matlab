function status = IntMetCom_Heterodyne(h,detector,ThreadID,BeamAID,BeamBID,...
    opdID, powID)
% status = IntMetCom_Heterodyne(...
%     h,detector,ThreadID,BeamAID,BeamBID,opdID,powID)
%
% heterodynes the beams stored in BeamAID and BeamBID
% detector.shape = 'circle';
% detector.size = [diameter or xsize, 0 or ysize]
% detector.xc, detector.yc
%
% results are stored in the COM Server outputOPD(opdID) and
% outputPOW(powID)
%
% at this time, there is no single-threaded call for heterodyne

% initialize to no error
status = 0;

rp = [detector.size(:)' detector.xc detector.yc];

CheckStatus(h.ClearThread(ThreadID));

CheckStatus(h.AddThreadCommand(ThreadID,IntMetComCMD('Heterodyne'),...
    rp,int32([opdID powID]),BeamAID,BeamBID,detector.shape));

status = h.ExecuteThread(ThreadID);

%%%%%%%%%%%%% nested function to check for error
    function CheckStatus(stmp)
        if stmp ~= 0,
            warning(['Heterodyne error: ' num2str(stmp)]);
            status = stmp;
            return
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end