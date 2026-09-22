function sOut = GenerateEFCReport_falco(runnum, TrialNum, listItnum, varargin)
% sOut = GenerateEFCReport_falco(runnum, TrialNum, listItnum, varargin)
%
% Backward-compatible wrapper around the CEFCReport class -- kept for
% existing callers (e.g. OMC_MSWC/run_testbed/GenReport.m). New code
% should use CEFCReport directly.
%
% listItnum can be array of itnums, or array of CRunData-subclass objects
% (e.g. CfalcoRunData)
% csDisplayFun is a cell array: {method, varargin options (e.g. 'clim', clim)}
% varargin = additional cdDisplayFun cell arrays, or bare flag names, as
% you want
%
% methods:
% DisplayImCubeImage
% DisplayImCubeUnProb
% DisplayImCubeContrast
% DisplayImCubeSigProb
% DisplayIncInt
% DisplayProbeAmp
% DisplayProbeCube
% DisplayCohInt
% DisplayAllInt
% DisplayRadialIntensity
% DisplayIncCohInt
% DisplayEfields
% DisplayDEfields
% DisplayCEfields
% DisplayDMv
% DisplayDMvProbe
%
% Options:
%   'pptfn', ''
%   'Sppt', 'new'            % 'new' | [] | existing Cppt object
%   'run_bn', 'falco_testbed_run<runnum>'
%   'listS', []              % existing array of CRunData-subclass, if iterations already loaded
%   'RunDataClass', 'CfalcoRunData'
%   'MaxIterGuess', 1000
%   'excludeItnum', []
%   ... plus display-plot flags/cells, see CEFCReport.AddDisplayPlots
%
% See also: CEFCReport

more off

listSin = CheckOption('listS', [], varargin{:}); % if listS of iterations already exists

% normalize listItnum: numeric itnums, or an array of CRunData-subclass objects
if isa(listItnum, 'CRunData')
    itnumRange = [listItnum.iter];
else
    itnumRange = listItnum;
end

if isempty(listSin)
    % "first call": fresh auto-discovery. The constructor builds the
    % default summary/compare/archived-figures sequence, and (since
    % itnumRange is forwarded as CEFCReport's own leading varargin
    % element) its internal AddDisplayPlots step renders exactly the
    % requested display plots restricted to itnumRange.
    R = CEFCReport(runnum, TrialNum, itnumRange, varargin{:});
else
    % "later call": reuse the already-loaded iteration data / open Sppt.
    % Do not re-run the summary/compare/archived-figures steps -- only
    % render the requested display plots.
    R = CEFCReport(runnum, TrialNum, varargin{:});
    R.AddDisplayPlots(itnumRange, varargin{:});
end

if nargout >= 1
    sOut = struct(...
        'listS', R.S ...
        ,'listHfig', {R.listHfig} ... % how to put a cell array in a struct field
        ,'Sppt', R.Sppt ...
        ,'probeh', {R.probeh} ...
        ,'rmsdDMv', {R.rmsdDMv} ...
        ,'itnum_best', R.itnum_best ...
        ,'fPlotNormIntensity', @R.PlotNormIntensity ...
        ,'fPlotBeta', @R.PlotBeta ...
        ,'fPlotProbeh', @R.PlotProbeh ...
        ,'fPlotRMSdDMv', @R.PlotRMSdDMv ...
        ,'fPlotCompareNormInt', @R.PlotCompareNormInt ...
        ,'Report', R ... % the CEFCReport object itself, for callers who want the OO style going forward
        );
end

more on

end % main
