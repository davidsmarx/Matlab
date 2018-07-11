function [n, dndt] = dispersionformula2index(p,formula,lam,d,T,Tref)
% [n, dndt] = dispersionformula2index(p,formula,lam,d,T,Tref)
%
% p is a row vector of coefficients, format depends on formula,
%   following the ZEMAX lens catalog format, see the ZEMAX manual
% formula is a string referring to one of the ZEMAX dispersion formulae.
% Valid choices are:
%   Sellmeier1: p = [k1 L1 k2 L2 k3 L3]
%            n^2 - 1 = k1 lam^2/(lam^2-L1) + k2 lam^2/(lam^2-L2) + ...
%   Schott: n^2 = p(1) + p(2)lam^2 + p(3)lam^-2 + p(4)lam^-4 +...
%   Herberger
%   Conrady: n = n0 + A/lam + B/lam^3.5
%   Cauchy: n = a + b/lam^2 + c/lam^4
%
% lam is wavelength (m)
% d (optional), coefficients for calculating dn/dt. Following the ZEMAX
% convention: d(1:3) = [D0 D1 D2], d(4:5) = [E0 E1], d(6) = lam_tk
% T (optional), temperature difference relative to reference temp.
% Tref (optional) usually = 20 C
%
% output n is relative to air (nabs = nair*nrel)
% output dndt is relative to air, use 'dndtschottformula' for abs dn/dT

global UM;

if nargin == 0, error('usage: n = dispersionformula2index(p,formula,lam,d,T)'); end

if isempty(formula), FormulaError; return; end

lam = lam(:)/UM; % formulae are all for lam in [um]

switch lower(formula)
    case 'schott'
        pp = [p(2) p(1) p(3:6)];
        n2 = (lam.^2).^-4.*polyval(pp,lam.^2);
        n = sqrt(n2);
        %keyboard;
        
    case 'sellmeier1'
        k = p([1:2:end]);
        l = p([2:2:end]);
        [K,Lam] = meshgrid(k,lam);
        [L,Lam] = meshgrid(l,lam);
        n2 = 1 + lam(:).^2.*sum(K./(Lam.^2 - L),2);
        n = sqrt(n2);
        
    case 'herzberger'
        L = 1./(lam(:).^2 - 0.028);
        n = (ones(length(lam(:)),1)*p(1:6)).*[ones(size(lam(:))) L L.^2 lam(:).^2 lam(:).^4 lam(:).^6];
        n = sum(n,2);        
        
    case 'conrady'
        n0 = p(1);
        A  = p(2);
        B  = p(3);
        n  = n0 + A./lam + B./(lam.^3.5);
        
    case 'cauchy'
        % n = a + b/lam^2 + c/lam^4 + ...
        lam = lam(:);
        n = repmat(1./lam,1,length(p)).^(2*repmat((0:length(p)-1),length(lam),1)) * p(:);
        
    otherwise,
        FormulaError;
        return
end


% if temperature is included, calculate index using temperature
if ~exist('T'), return, end
if ~exist('Tref'), Tref = 20; end

T = T(:);

%Tref = d(7);
dT = T - Tref;

nairref = indexofair(lam,Tref,1)';
nabsref = n .* nairref;

dn = d(1)*dT + d(2)*dT.^2 + d(3)*dT.^3 + (d(4)*dT + d(5)*dT.^2)./(lam.^2 - d(6)^2);
dn = 0.5*(n.^2-1).*dn./n;

nabs = nabsref + dn;

[nair, dndtair] = indexofair(lam,T,1);

% if requested, calculate dn/dt
if nargout >= 2,
   dndtabs = dndtschottformula(d,n,lam,dT(:));
   dndt = (1./nair).*(dndtabs - n.*dndtair);
end

n = nabs(:)./nair(:);

    function FormulaError
        disp(['formula ' formula ' is not valid']);
        n = NaN;
        dndt = NaN;
    end

end
