function [Hb, Ha, Gb, Ga] = rcoeffs2tf(rf)
% [Hb, Ha, Gb, Ga] = rcoeffs2tf(rf,rb,tf,tb)
% rf = reflection coefficients at each thin film boundary
% r(1) = substrate-to-first layer boundary
% Hb = numerator polynomial of z-transform of reflection tf
% Ha = denominator polynomial of z-transform of reflection tf
% gb = numerator polynomial of z-transform of transmission tf
% ga = denominator
if nargin == 0, error('usage: [Hb, Ha, Gb, Ga] = rcoeffs2tf(rf,rb,tf,tb)'); end

Nr = length(rf);

syms z;

PQ = [rf(1); 1];
for i = 2:Nr,
   PQ = [1 rf(i)*z; conj(rf(i)) z] * PQ;

end

Hb = PQ(1);
Ha = PQ(2);

Gb = prod(sqrt(1 - abs(rf).^2));
Ga = PQ(2);

return