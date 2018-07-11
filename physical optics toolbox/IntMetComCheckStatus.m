function IntMetComCheckStatus(status)
% IntMetComCheckStatus(status)

% INTEGER, PARAMETER :: INVALIDTHREADID = -10
% INTEGER, PARAMETER :: INVALIDWFRONTID = -11
% INTEGER, PARAMETER :: INVALIDOUTPUTID = -12
% !! constant definitions for flagStatus
% INTEGER, PARAMETER :: STATUS_NOERROR    =  0
% INTEGER, PARAMETER :: STATUS_WORKING    =  1
% INTEGER, PARAMETER :: STATUS_EXITTHREAD = -1
% INTEGER, PARAMETER :: STATUS_ERROR      = -2
% INTEGER, PARAMETER :: STATUS_CMDIDERROR = -3

switch status,
    case 0,
        % no error, do nothing
    case 1,
        % STATUS_WORKING, do nothing
    case -1,
        error('STATUS_EXITTHREAD');
    case -2,
        error('STATUS_ERROR');
    case -3,
        error('STATUS_CMDIDERROR');
    case -10,
        error('INVALIDTHREADID');
    case -11,
        error('INVALIDWFRONTID');
    case -12,
        error('INVALIDOUTPUTID');
    otherwise,
        error(['unknown error: ' num2str(status)]);
end
