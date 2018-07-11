function retval = z_setwave(zchan,lam,nwave,weight)
% status = z_setwave(zchan,lam,nwave,weight)
% lam = wavelength (m)
% nwave (optional) = integer wavelength number (default = 1)
% weight (optional) = weight applied to this wavelength (default = 1.0)
%
% status: if nwave ~= 0, status = [wavelength weight]
%         if nwave == 0, status = [primary wavelength #,
%                                  number of wavelengths currently defined]
%
% GetUpdate is called after setting the new wavelength

global UM;

if nargin == 0, disp('usage: status = z_setwave(zchan,lambda)'); return, end

if ~exist('nwave','var') | isempty(nwave),
    nwave = 1;
end
if ~exist('weight','var') | isempty(weight),
    weight = 1.0;
end

% ZEMAX uses um for wavelength
cmdstr = sprintf('SetWave,%d,%f,%f',nwave,lam/UM,weight);
stmp = ddereq(zchan,cmdstr,[1 1]);
if ~ischar(stmp), error(['SetWave error, return value = ' num2str(stmp)]); end

retval = sscanf(stmp,'%f,',[1 inf]);

% call update
%if z_getupdate(zchan), error('z_setwave: GetUpdate Failed!'); end

return
