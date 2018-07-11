function [ retval ] = z_NSCTrace(zchan, sParms )
%  retval = z_NSCTrace(zchan, sParms )
%
% Clear Detectors, and then trace
%
% NSCTrace,surf,source,split,scatter,usepolar,ignore_errors,no_random_seed,save,
% savefilename,filter,zrd_format
%
% sParms = struct(...
%     'surf', 1 ... % if mode is set to Non-sequential, use 1
%     ,'source', 0 ... % trace all sources
%     ,'split', 1 ... % non-zero => splitting is on
%     ,'scatter', 0 ...% 0 => scattering is off
%     ,'usepolar', 1 ... % use polarization
%     ,'ignore_errors', 0 ... % errors will terminate the ray trace
%     ,'random_seed', 0 ... % 0 => each ray trace uses a different random seed, when using NSTR for optimization, it is recommended to use non-zero seed
%     ,'save', 0 ...% don't save the ray data
%     ,'savefilename', '""' ...
%     ,'filter', '""' ...
%     ,'zrd_format', '""' ...
%     );

% check parameters
if nargin == 0,
    % just return the default sParms struct
    retval = ValidateParms(struct());
    return
end

if ~exist('sParms','var'), sParms = struct; end
sParms = ValidateParms(sParms);

% be sure to clear detectors using:  z_getNSCDetectorData(ZEMAXDDE);
% % clear detectors
% % NSDD 1, 0 (clear all detectors), 0, 0
% % NSCDetectorData,surface,object,pixel,data
% cmdstr = 'NSCDetectorData,1,0,0,0';
% stmp = ddereq(zchan,cmdstr,[1 1]);


% make the call
cmdstr =...
    ['NSCTrace' ...
    ', ' num2str(sParms.surf) ...
    ', ' num2str(sParms.source) ...
    ', ' num2str(sParms.split) ...
    ', ' num2str(sParms.scatter) ...
    ', ' num2str(sParms.usepolar) ...
    ', ' num2str(sParms.ignore_errors) ...
    ', ' num2str(sParms.random_seed) ...
    ', ' num2str(sParms.save) ...
    ', ' sParms.savefilename ...
    ', ' sParms.filter ...
    ', ' sParms.zrd_format ...
    ];

if length(cmdstr) > 255,
    error(['length of cmdstr = ' num2str(length(cmdstr))]);
end

    retval = ddereq( zchan,cmdstr,[1 1]);

% wait until zemax is ready
while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end

end

function sParms = ValidateParms(sPin)

% default values: see p. 770 of manual
sParms = struct(...
    'surf', 1 ... % if mode is set to Non-sequential, use 1
    ,'source', 0 ... % trace all sources
    ,'split', 1 ... % non-zero => splitting is on
    ,'scatter', 0 ...% 0 => scattering is off
    ,'usepolar', 1 ... % use polarization
    ,'ignore_errors', 0 ... % errors will terminate the ray trace
    ,'random_seed', 0 ... % 0 => each ray trace uses a different random seed, when using NSTR for optimization, it is recommended to use non-zero seed
    ,'save', 0 ...% don't save the ray data
    ,'savefilename', '""' ...
    ,'filter', '""' ...
    ,'zrd_format', '""' ...
    );

% copy input parameters
fnames = fieldnames(sPin);
for ii = 1:length(fnames),
    sParms.(fnames{ii}) = sPin.(fnames{ii});
end

end