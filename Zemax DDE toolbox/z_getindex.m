function [nn, lam] = z_getindex(zchan,ss)
% nn = z_getindex(zchan,ss)
% ss = either surface struct, or array of surface numbers
% nn = array of indexes, 
%      one column for each surface, one row for each wavelength
% lam = (optional output argument) list of defined wavelengths

if nargin == 0, error('usage: [nn, lam] = z_getindex(zchan,ss)'); end

if isstruct(ss), nsurf = [ss.nsurf]; else, nsurf = ss(:)'; end

nn = [];
for ns = nsurf,
   cmdstr = ['GetIndex, ' num2str(ns)];
   a = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);
   nn = [nn a(:)];
end

if nargout > 1,
   for i = 1:length(a),
      cmdstr = ['GetWave, ' num2str(i)];
      b = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);
      lam(i) = b(1);
   end
end


return


