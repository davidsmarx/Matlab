function sOut = GenerateEFCReport(runnum, listItnum, sDisplayFun, varargin)
% sOut = GenerateEFCReport(runnum, listItnum, sDisplayFun)
%
% listItnum can be array of itnums, or array of CRunData objects
%
% create figures S.(sDisplayFun) for each iteration in listItnum
% and copy to a PowerPoint presentation
mlock
   

persistent Sppt;

% options
ppt_fn = CheckOption('pptfn', '', varargin{:});
bDiff  = CheckOption('differential', false, varargin{:}); %

% open PowerPoint if necessary
if isempty(Sppt)
    Sppt = Cppt(ppt_fn);
end

% get the CRunData objects
N = length(listItnum);
if isnumeric(listItnum),
    for ii = 1:N
        fprintf('reading itnum %d\n', listItnum(ii));
        S(ii) = CRunData(runnum, listItnum(ii));
    end
elseif isa(listItnum, 'CRunData')
    S = listItnum;
else
    error(['listItnum type error: ' class(listItnum)]);
end
  
% create the plots
% some plots are differential
if ~bDiff,
for ii = 1:N,    
    hfig(ii) = S(ii).(sDisplayFun)(varargin{:});
    set(hfig(ii), 'Position', 0.6*get(hfig(ii),'position'));
    Sppt.CopyFigNewSlide(hfig(ii));
end

else
for ii = 1:N-1,
    hfig(ii) = S(ii+1).(sDisplayFun)(S(ii), varargin{:});
    set(hfig(ii), 'Position', 0.6*get(hfig(ii),'position'));
    Sppt.CopyFigNewSlide(hfig(ii));
end    
end

if nargout >= 1,
    sOut = struct(...
        'listS', S ...
        ,'listHfig', hfig ...
        );
end