function tfperf = TFfilter_performance(lam,il,options)
% tfperf = TFfilter_performance(lam,il,options)
%
% return struct:
%   tfperf.cfr     [THz] center frequency scalar
%   tfperf.frbw    [THz] bandwidth at given level, vector size of levdb
%   tfperf.lambw   [um]  same as frbw
%   tfperf.pbil    [dB]  maximum insertion loss over the pass band (+/-bwfr)
%   tfperf.pbrip
%   tfperf.leftadjrej
%   tfperf.rigtadjrej
%   tfperf.lnonadjrej
%   tfperf.rnonadjrej
%   tfperf.fr     [THz] freq grid to calculate gd, cd, not output if input lam=freq.
%   tfperf.phi    [rad]
%   tfperf.gd     [ps] sampled on the frequency grid (not wavelength)
%   tfperf.cd     [ps/nm] sampled on the frequency grid
%   tfperf.gdrip  [ps], max(GD)-min(GD) within +/-bwfr
%   tfperf.cdmax  [ps/nm], max |CD| within +/-bwfr
% input:
%   lam = wavelength [um] (default) or frequency [THz]
%   il [dB] (default) or [linear complex amplitude]
%      if il is complex, then GD and CD are calculated directly,
%      otherwise, GD and CD are calculated from minimum phase response
% options is a struct with the following optional fields:
%   'freq' : input lam is frequency [THz] (default is wavelength)
%   'reflected' : input il is for reflected channel (default is transmitted)
%   'scale' : 'dB', input il is [dB] (default)
%             'linear', input il is linear power scale
%             'complex', input il is linear complex, GD and CD will be directly calculated
%   'gd' : vector of measured GD (optional), CD is then calculated from the gradient of GD.
%          scale must be either 'linear' or 'dB'. If gd data is not provided,
%          then gd is calculated from the minimum phase response.
%   'levdb' followed by scalar [dB] vector of levels at which to calculate
%         bandwidth (bw), cwl is calculated at levdb(1)
%   'bwfr' followed by a scalar [THz] : half-bandwidth (bandwidth is +/- bwfr)
%         used for IL, ripple, and max CD
%   'adjfr' followed by a scalar [THz] : frequency offset to edge of adjacent channel
%   'nonadjfr' followed by a scalar [THz] : frequency offset to non adjacent channel
%   'cfr' : specified center frequency (e.g. ITU grid) [THz] (required for reflected)
%   'chsepfr' : channel separation [THz] (required for reflected)

if nargin == 0,
   error('usage: tfperf = TFfilter_performance(lam,il,options)');
end
if ~exist('options','var'), options.scale = 'dB'; end

constants;

if isfield(options,'freq'),
   fr = lam(:); % input is frequency
else, % convert wavelength to frequency and resample il
   fr  = linspace(C./lam(end),C./lam(1),length(lam))'; % [THz]
   il  = pchip(C./lam,il,fr);
   tfperf.fr = fr;
end

[tfperf.gd, tfperf.cd, ildb] = calculate_gdcd(lam,fr,il,options);

if isfield(options,'reflected'),
   %%%%%%%%%%% reflected response performance %%%%%%%%%%%%%%
   if ~isfield(options,'cfr') | ~isfield(options,'chsepfr'),
      error('center frequency of passband (cfr) and channel spacing (chsepfr) must be defined');
   end
   if isfield(options,'bwfr'),
      [tfperf.chrej, tfperf.ril, tfperf.cdleftadj, tfperf.cdrigtadj] =...
         calc_reflparms(fr,ildb,tfperf.cd,options);
   end
   
else,
   %%%%%%%%%%% transmitted response performance %%%%%%%%%%%%
   [tfperf.cfr, tfperf.lambw, tfperf.frbw] = calculate_bw(fr,ildb,options);

   % use required center frequency for passband and adjacent band performance
   if isfield(options,'cfr'),
      cfr = options.cfr;
   else,
      cfr = tfperf.cfr;
   end
   
   if isfield(options,'bwfr'),
      [tfperf.pbil, tfperf.pbrip, tfperf.gdrip, tfperf.cdmax] =...
         calc_rips(fr,cfr,ildb,tfperf.gd,tfperf.cd,options);
   end
   
   % adjacent channel rejection
   if isfield(options,'adjfr'),
      [tfperf.leftadjrej, tfperf.rigtadjrej] = calc_adj(fr,cfr,ildb,options);
   end
   
   % non-adjacent channel rejection
   if isfield(options,'nonadjfr'),
      nonadjfr = options.nonadjfr;
      tfperf.lnonadjrej = interp1(fr,ildb,cfr(1)-nonadjfr) - max(il);
      tfperf.rnonadjrej = interp1(fr,ildb,cfr(1)+nonadjfr) - max(il);
   end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       routines for calculating various performance parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gd,cd,ildb] = calculate_gdcd(lam,fr,il,options)

% if il represents complex amplitude signal, then convert to dB
% calculate GD and CD
constants;
if ~isfield(options,'scale'), options.scale = 'dB'; end
switch options.scale,
case 'complex',
   % direct calculation of GD and CD from phase
   [gd, cd] = group_delay(il,diff(fr([1 2])),mean(fr));
   tfperf.phi = angle(il);
   ildb = 2*decibel(il);
   
case {'linear', 'dB'},
   if strcmp(options.scale,'linear'),
      ildb = decibel(il);
   else,
      ildb = il;
   end
   
   if isfield(options,'gd'),
      % measured data includes IL and GD (input must be wavelength)
      % cd is [ps/nm] if lam is [um] and gd is [ps]
      
      % resample gd onto even wavelength grid to take gradient   
      lami = linspace(C./fr(end),C./fr(1),length(fr));
      if isfield(options,'freq'), % gd is function of input freq (lam)
         gd = pchip(C./lam,options.gd,lami);
      else, % gd is function of input lam
         gd = pchip(lam,options.gd,lami);
      end

      % calculate gradient on wavelength grid, then return to frequency grid
      % filter measurement data
%      b = remez(28,[0 0.05 0.1 1],[1 1 0 0]); b = b./sum(b);
%       b = remez(28,[0 0.2 0.25 1],[1 1 0 0]); b = b./sum(b);
      b = remez(28,[0 0.4 0.5 1],[1 1 0 0]); b = b./sum(b);
      gdfilt = filtfilt(b,1,gd);
      cd = gradient(gdfilt,1000*diff(lami([1 2])));
%       cd = gradient(gd,1000*diff(lami([1 2])));

      gd = pchip(C./lami,gdfilt,fr);
%       gd = pchip(C./lami,gd,fr);
      cd = pchip(C./lami,cd,fr);
      
   else,   
      % calculate min phase for measured IL, passband CD
      [phi, gd, cd] = min_phase(fr,il,'THz','dB','um','Ngrad','16','Fgrad','0.1');
      tfperf.phi = phi;
   end
      
otherwise,
   error('unknown scale type');
end

return

function [cfr, lambw, frbw] = calculate_bw(fr,il,options)
% calculate filter bandwidth for given level (e.g. -1dB bandwidth)
constants;
if isfield(options,'levdb'),
   levdb = options.levdb;
   ilmax = max(il);
   for ii = 1:length(levdb),
      [cfr(ii),frbw(ii)] = cwlbw(fr,il-ilmax,levdb(ii),'dB');
      lambw(ii) = C./(cfr(ii)-0.5*frbw(ii)) - C./(cfr(ii)+0.5*frbw(ii));
      % correct bw for cfr shift, use levdb(1) as reference center
      frbw(ii) = frbw(ii) + abs(cfr(ii)-cfr(1));
   end
%    frbw = frbw(:);
%    lambw = lambw(:);
else,
   [ilmax, imax] = max(il);
   cfr = fr(imax);
   lambw = [];
   frbw = [];
end
cfr = cfr(1);

return

function [pbil, pbrip, gdrip, cdmax] = calc_rips(fr,cfr,il,gd,cd,options)
% parameters based on given frequency passband (e.g. ripple +/-10GHz)

for ii = 1:length(options.bwfr),
   bwfr = options.bwfr(ii);
   
   npass = fr > cfr(1)-bwfr & fr < cfr(1)+bwfr;
   % IL, ripple, group delay ripple, and max CD
   pbil(ii) = max(abs(il(npass)));
   pbrip(ii) = max(il(npass)) - min(il(npass));
   gdrip(ii) = max(gd(npass)) - min(gd(npass));
   
   cdpass = cd(npass);
   [tmpmax, imax] = max(abs(cdpass));
   cdmax(ii) = cdpass(imax);
end

return

function [leftadjrej, rigtadjrej] = calc_adj(fr,cfr,il,options)
   
   adjfr = options.adjfr;
   nleftadj = fr < cfr-adjfr;
   nrigtadj = fr > cfr+adjfr;
   leftedge = interp1(fr,il,cfr-adjfr) - max(il);
   rigtedge = interp1(fr,il,cfr+adjfr) - max(il);
   leftadjrej = max([leftedge; il(nleftadj)]);
   rigtadjrej = max([rigtedge; il(nrigtadj)]);

return

function [chrej, ril, cdleftadj, cdrigtadj] = calc_reflparms(fr,il,cd,options);
% parameters for reflected channel rejection, RIL, and adjacent channel cd

npass = fr > options.cfr - options.bwfr & fr < options.cfr + options.bwfr;
chrej = max(il(npass));

nleftrefl = fr < options.cfr - options.chsepfr + options.bwfr;
nrigtrefl = fr > options.cfr + options.chsepfr - options.bwfr;
ril = min([min(il(nleftrefl)) min(il(nrigtrefl))]);

cdleft = cd(nleftrefl);
[cdmax, imax] = max(abs(cdleft));
cdleftadj = cdleft(imax);

cdrigt = cd(nrigtrefl);
[cdmax, imax] = max(abs(cdrigt));
cdrigtadj = cdrigt(imax);

return
