function [TDone, TNotDone] = IntMetCom_WaitForThreads(h,Tlist,TimeOut,flagWaitAll)
% [TDone, TNotDone] = IntMetCom_WaitForThreads(h,Tlist,TimeOut,flagWaitAll)
%
% h = COM server
% Tlist = vector of Thread ID's
% TimeOut = seconds (default = 60)
%
% flagWaitAll (optional): true (default) = wait for all threads to finish,
% false = return when any one thread is finished.
%
% TDone = list of threadID's that signaled
% TNotDone = list of threadID's that are not signaled

global MS;

% check for optional TimeOut
if ~exist('TimeOut','var') | isempty(TimeOut),
    TimeOut = 60;
end

% check for optional flagWaitAll
if ~exist('flagWaitAll','var') | isempty(flagWaitAll),
    flagWaitAll = 1;
end

nT = length(Tlist);

if nT == 1,
    waitstatus = h.WaitForThread(Tlist,TimeOut/MS);
elseif nT > 1,
    waitstatus = h.WaitForMultipleThreads(int32(Tlist),TimeOut/MS,flagWaitAll);
else, % must be empty
    error('empty thread list, Tlist');
end
CheckWaitStatus;

%%%%%%%%% from DFWINTY.F90
%   integer, parameter :: WAIT_FAILED = ( #FFFFFFFF)
%   integer, parameter :: WAIT_OBJECT_0 = ((( #00000000) ) + 0 )
%   integer, parameter :: WAIT_ABANDONED = ((( #00000080) ) + 0 )
%   integer, parameter :: WAIT_ABANDONED_0 = ((( #00000080) ) + 0 )
%   integer, parameter :: WAIT_TIMEOUT = 258

    function CheckWaitStatus
        switch waitstatus,
            case mat2cell([0:nT-1],1,ones(1,nT)),
                % thread terminated successfully (signaled state)
                if flagWaitAll,
                    TDone = Tlist;
                    TNotDone = [];
                else,
                    iT = waitstatus + 1;
                    TDone = Tlist(iT);
                    TNotDone = Tlist([1:iT-1 iT+1:end]);
                end
                CheckError(TDone);
            case 258,
                % timeout
                warning(['Time Out Error']);
                keyboard;
                TDone = [];
                TNotDone = Tlist;
            case mat2cell(128+[0:nT-1],1,ones(1,nT)),
                % wait abandoned
                error(['Wait Abandoned, ThreadID = ' num2str(Tlist(waitstatus-127))]);
            case -1,
                % wait function failed, call GetLastError
                error(['Wait Failed ' datestr(now)]);
            otherwise,
                error(['Unknown wait status = ' num2str(waitstatus)]);
        end
        
    end
%%%%%%%%%%%%%%%%%%%%
    function CheckError(Thread)
        for ii = 1:length(Thread),
            if h.ThreadStatus(Thread(ii)) < 0,
                error(['thread error: ' num2str(h.ThreadStatus(Thread(ii)))]);
            end
        end
    end
%%%%%%%%%%%%%%%%%%%
end