function tfperf = TFfilter_performance(lam,il,varargin)
% tfperf = TFfilter_performance(lam,il,levdb,bwfr,adjfr,nonadjfr)
%levdb,bwfr,adjfr,nonadjfr)
% return struct:
%   tfperf.cfr     [THz] center frequency scalar
%   tfperf.frbw    [THz] bandwidth at given level, vector size of levdb
%   tfperf.lambw   [um]  same as frbw
%   tfperf.pbrip
%   tfperf.leftadjrej
%   tfperf.rigtadjrej
%   tfperf.lnonadjrej
%   tfperf.rnonadjrej
%   tfperf.phi    [rad]
%   tfperf.gd     [ps]
%   tfperf.cd     [ps/nm]
%   tfperf.gdrip  [ps], max |GD| within +/-bwfr
%   tfperf.cdmax  [ps/nm], max |CD| within +/-bwfr
% input:
%   lam = wavelength [um] (default) or frequency [THz]
%   il [dB] (default) or [linear complex amplitude]
%      if il is complex, then GD and CD are calculated directly,
%      otherwise, GD and CD are calculated from minimum phase response
% options:
%   'freq' : input lam is frequency [THz] (default is wavelength)
%   'dB' : input il is [dB] (default)
%   'linear' : input il is linear power scale
%   'complex' : input il is linear complex, GD and CD will be directly calculated
%   'levdb' followed by scalar [dB] vector of levels at which to calculate
%         bandwidth (bw), cwl is calculated at levdb(1)
%   'bwfr' followed by a scalar [THz] : half-bandwidth (bandwidth is +/- bwfr)
%         used for ripple and max CD
%   'adjfr' followed by a scalar [THz] : frequency offset to edge of adjacent channel
%   'nonadjfr' followed by a scalar [THz] : frequency offset to non adjacent channel

if nargin == 0,
   error('usage: tfperf = TFfilter_performance(lam,il,levdb,bwfr,adjfr,nonadjfr)');
end
constants;

if any(strcmp('freq',varargin)),
   fr = lam(:); % input is frequency
else, % convert wavelength to frequency and resample il
   fr  = linspace(C./lam(end),C./lam(1),length(lam))'; % [THz]
   il  = pchip(C./lam,il,fr);
end

% if il represents complex amplitude signal, then convert to dB
% calculate GD and CD
if any(strcmp('complex',varargin)),
   % direct calculation of GD and CD from phase
   [gd, cd] = group_delay(il,diff(fr([1 2])),mean(fr));
   tfperf.phi = angle(il);
   il = 2*decibel(il);
   
else,
   if any(strcmp('linear',varargin)),
      il = decibel(il);
   end
   % calculate min phase for measured IL, passband CD
   [phi, gd, cd] = min_phase(fr,il,'THz','dB','um','Ngrad','16','Fgrad','0.1');
   tfperf.phi = phi;
end
tfperf.gd = gd;
tfperf.cd = cd;
   
% calculate filter bandwidth for given level (e.g. -1dB bandwidth)
if any(strcmp('levdb',varargin)),
   levdb = varargin{find(strcmp('levdb',varargin))+1};
   ilmax = max(il);
   for ii = 1:length(levdb),
      [cfr(ii),frbw(ii)] = cwlbw(fr,il-ilmax,levdb(ii),'dB');
      lambw(ii) = C./(cfr(ii)-0.5*frbw(ii)) - C./(cfr(ii)+0.5*frbw(ii));
   end
   tfperf.frbw = frbw(:);
   tfperf.lambw = lambw(:);
else,
   [ilmax, imax] = max(il);
   cfr = fr(imax);
end
tfperf.cfr = cfr(1);

% parameters based on given frequency passband (e.g. ripple +/-10GHz)
if any(strcmp('bwfr',varargin)),
   bwfr = varargin{find(strcmp('bwfr',varargin))+1};
   npass = fr > cfr(1)-bwfr & fr < cfr(1)+bwfr;
   % ripple, group delay ripple, and max CD
   tfperf.pbrip = max(il(npass)) - min(il(npass));
   tfperf.gdrip = max(gd(npass)) - min(gd(npass));
   cdpass = cd(npass);
   [cdmax, imax] = max(abs(cdpass));
   tfperf.cdmax = cdpass(imax);
end

% adjacent channel rejection
if any(strcmp('adjfr',varargin)),
   adjfr = varargin{find(strcmp('adjfr',varargin))+1};
   nleftadj = fr < cfr(1)-adjfr;
   nrigtadj = fr > cfr(1)+adjfr;
   leftedge = interp1(fr,il,cfr(1)-adjfr) - ilmax;
   rigtedge = interp1(fr,il,cfr(1)+adjfr) - ilmax;
   tfperf.leftadjrej = max([leftedge; il(nleftadj)]);
   tfperf.rigtadjrej = max([rigtedge; il(nrigtadj)]);

end

% non-adjacent channel rejection
if any(strcmp('nonadjfr',varargin)),
   nonadjfr = varargin{find(strcmp('nonadjfr',varargin))+1};
   tfperf.lnonadjrej = interp1(fr,il,cfr(1)-nonadjfr) - ilmax;
   tfperf.rnonadjrej = interp1(fr,il,cfr(1)+nonadjfr) - ilmax;
end

return
