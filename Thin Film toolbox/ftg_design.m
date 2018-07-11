function retstr = ftg_design(designstr)
% retstr = ftg_design(designstr)
% designstr example: '0.25H 0.25L 0.25H'
% if designstr = [], then just return the current design
% return design if successful, [], if not

% First open Server
Design = actxserver('Ftgdesign1.clsBasic');

try % if anything goes wrong, be sure to release Design
	
	% 
	if ~isempty(designstr),
		clipboard('copy',designstr);
		invoke(Design,'DesignPaste');
	end
	
	retstr = Design.Design;
	
	% execute spectrum calculation
	invoke(Design,'CalcPlot');
	
	% get spectrum
	lam = Design.Spectrum_X;
	spectrum = num2cell(double(Design.Spectrum_Y),1);

catch,
	retstr = [];
end
release(Design);
	
switch length(spectrum),
case 2, 
	[r, t] = deal(sqrt(spectrum{:})); % intensity only
	
case 4,
	[Rphi, R, Tphi, T] = deal(spectrum{:}); % phase and intensity
	r = sqrt(R).*exp(-j*Rphi);
	t = sqrt(T).*exp(-j*Tphi);
otherwise,
	error('unkown spectrum data type');
end

plotampphase(lam,[r t],'dB','legend','Reflected','Transmitted');

% convert lam back to [um] before returning
lam = 0.001*double(lam);

return