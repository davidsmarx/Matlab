function [r, t, lam] = ftg_spectrum(lam1,lam2,N,varargin)
% [r, t, lam] = ftg_spectrum(lam1,lam2,N,[optionlist])
% lam1, lam2 are in [um]
% lam1, lam2 can be aoi in [deg]
% supported options:
%    'tffdesign' string containing layer thickness, substrate first, in units of 
%             quarter wave. e.g. 1L 1H 1L
%    'pol' 'R'andom, 'P' or 'S'
%    'aoi' 'degrees as a string'
%    'plotaoi' calculate response vs. aoi, set lam1 and lam2 to aoi [deg]
%    'evalwave' evaluation wavelength as a string [nm], for use when 'plotaoi'
%    'designwave' design wavelength [nm] as string
%    'optimize' causes Film Star to run optimize before calculating spectrum
%
% return complex field coefficients and a vector of wavelengths [um]


% First open Server
Design = actxserver('Ftgdesign1.clsBasic');

try % if anything goes wrong, be sure to release Design

   ides = strmatch('tffdesign',varargin);
	if ides,
		clipboard('copy',varargin{ides+1});
		invoke(Design,'DesignPaste');
	end

   ipol = strmatch('pol',varargin);
	if ipol, 
		if any(varargin{ipol+1}(1) == 'PSR'),
			Design.Pol = varargin{ipol+1}(1);
		else,
			warning('invalid polarization setting, not changing polarization');
		end
	end
	
	iaoi = strmatch('aoi',varargin);
	if iaoi,	Design.Angle = varargin{iaoi+1};	end
	
	if strmatch('plotaoi',varargin),
		flam = 1;
	else,
		% convert lam to [nm]
		flam = 1000; lam1 = flam*lam1; lam2 = flam*lam2;
	end

	iewa = strmatch('evalwave',varargin);
	if iewa, Design.EvalWave = varargin{iewa+1}; end
	
	idwa = strmatch('designwave',varargin);
	if idwa, Design.DesignWave = varargin{idwa+1}; end
	
	% setup Film Star Spetrum axes
	axesstr = sprintf('%f\t%f\t%f\t%d\n0\t30\t5\nNatick',...
		lam1,lam2,abs(lam2-lam1)/(N-1),ceil(abs(lam2-lam1)/5));
	clipboard('copy',axesstr);
	invoke(Design,'AxesPaste');

   % run optimization before calculating spectrum?
   iopt = strmatch('optimize',varargin);
   if iopt, invoke(Design,'optimize'); end
   
	% execute spectrum calculation
	invoke(Design,'Calculate');
	
	% get spectrum
	lam = double(Design.Spectrum_X)' ./ flam;
	spectrum = num2cell(double(Design.Spectrum_Y),1);
   
catch,
	warning('film star error');
	keyboard;
end
release(Design);
	
switch length(spectrum),
case 2, 
	[r, t] = deal(spectrum{:}); % intensity only
	r = sqrt(r);
	t = sqrt(t);
	
case 4,
	[Rphi, R, Tphi, T] = deal(spectrum{:}); % phase and intensity
	r = sqrt(R).*exp(-j*Rphi);
	t = sqrt(T).*exp(-j*Tphi);
otherwise,
	warning('unkown spectrum data type');
	keyboard;
end

%plotampphase(lam,[r t],'dB','legend','Reflected','Transmitted');

return