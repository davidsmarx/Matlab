function [xm1, am1, xm2, am2] = find2peaks(x, amp, parmstruct)
% [xm1, am1, xm2, am2] = find2peaks(x, amp, parameters
% 
% parameters struct fields:
%    minAmp (default = 0)
%    minWidth (default = 1)
%    minAlpha (defualt = 0.1)
%    FilterLength (default = 11)

% check input parameters
if exist('parmstruct','var'),
    if isstruct(parmstruct),
        parms = CheckParms(parmstruct);
    else,
        error('invalid parameters');
    end
else,
    parms = DefaultParms;
end

samplespacing = mean(diff(x));
NDKERN = floor(parms.FilterLength/samplespacing);
NDKERN_2 = floor(NDKERN/2);

% filter the input signal
% b = [ones(1,NDKERN)]/NDKERN;
% a = 1;
% amp = circshift(filter(b,a,amp),-NDKERN_2);
% figure, plot(x,amp), grid

% apply derivative filter
b = [ones(1,NDKERN_2) 0 -ones(1,NDKERN_2)];
a = 1;
derv = circshift(filter(b,a,amp),-NDKERN_2);

% zero crossings and peak locations
zcrossings = [0; diff(derv > 0)];
ipeaks = find(zcrossings < 0);
ivalleys = find(zcrossings > 0);

%disp(ipeaks);
%disp(ivalleys);
peaks = x(ipeaks);

% for each peak, get the closest valley to the right and left
% get peak amplitude
for ii = 1:length(ipeaks),
    ip = ipeaks(ii);
    ivltmp = max(ivalleys(ivalleys < ip));
    ivrtmp = min(ivalleys(ivalleys > ip));
    if isempty(ivltmp) || isempty(ivrtmp),
        peakamp(ii) = 0;
        peakwidth(ii) = 0;

    else,
        ivalleyleft(ii) = ivltmp;
        ivalleyright(ii) = ivrtmp;
        ampL(ii) = amp(ip) - amp(ivalleyleft(ii));
        ampR(ii) = amp(ip) - amp(ivalleyright(ii));
        peakamp(ii) = max(ampL(ii),ampR(ii));
        peakwidth(ii) = x(ivalleyright(ii)) - x(ivalleyleft(ii));

    end
    
end %for each peak

%[ipeaks(:) ivalleyleft(:) ivalleyright(:) ampL(:) ampR(:) peakamp(:) peakwidth(:)]

% pick the strongest peak
[maxamp, i1] = max(peakamp);
[am1, xm1] = findpeakn(x(ipeaks(i1)+[-NDKERN_2:NDKERN_2]),amp(ipeaks(i1)+[-NDKERN_2:NDKERN_2]),1);

% elliminate that peak and test the rest against criteria
peakamp(i1) = 0;
peakamp(peakamp < parms.minAmp) = 0;
peakamp(peakwidth < parms.minWidth) = 0;
peakamp(peakamp < parms.minAlpha*am1) = 0;

% now get second peak from what's left
[maxamp, i2] = max(peakamp);
if maxamp > 0,
    [am2, xm2] = findpeakn(x(ipeaks(i2)+[-NDKERN_2:NDKERN_2]),amp(ipeaks(i2)+[-NDKERN_2:NDKERN_2]),1);
else,
    xm2 = xm1;
    am2 = am1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parms = DefaultParms
parms.minAmp = 0;
parms.minWidth = 1;
parms.minAlpha = 0.1;
parms.FilterLength = 11;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parms = CheckParms(parmstruct)
defparms = DefaultParms;

parms = parmstruct;

fnames = fieldnames(defparms);
for ii = 1:length(fnames),
    if ~isfield(parmstruct,fnames{ii}),
        parms.(fnames{ii}) = defparms.(fnames{ii});
    end
end
end
