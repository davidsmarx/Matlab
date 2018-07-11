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

P = rf(1); Q = 1;
for i = 2:Nr,
   Pi = rf(i)*[Q 0] + [0 P];
   Qi = [Q 0] + conj(rf(i))*[0 P];
   
   P = Pi;
   Q = Qi;
end

Hb = P;
Ha = Q;

Gb = prod(sqrt(1 - abs(rf).^2));
Ga = Q;

return