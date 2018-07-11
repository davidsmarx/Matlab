function [Hb, Ha, Gb, Ga] = rcoeffs2tf(rf,rb,tf,tb)
% [Hb, Ha, Gb, Ga] = rcoeffs2tf(rf,rb,tf,tb)
% rf = reflection coefficients at each thin film boundary
%   optional:
%      rb = reflection coefficients at each thin film boundary
%           default: rb = rf
%      tf = forward transmission coeff. (incident med. to sub.)
%           default: tf = sqrt(1 - rf.^2)
%      tb = backward transmission coeff. (sub. to inc. med.)
%           default: tb = tf
% r(1), t(1) = substrate-to-first layer boundary
% Hb = numerator polynomial of z-transform of reflection tf
% Ha = denominator polynomial of z-transform of reflection tf
% gb = numerator polynomial of z-transform of transmission tf
% ga = denominator
if nargin == 0, error('usage: [Hb, Ha, Gb, Ga] = rcoeffs2tf(rf,rb,tf,tb)'); end

Nr = length(rf);

if nargin <= 2,
	tf = sqrt(1 - rf.^2);
	tb = tf;
   if nargin == 1, rb = rf; end
elseif nargin == 4,
	if length(rb) ~= Nr | length(tf) ~= Nr, error('r and t must be same size'); end
else,
	error('incorrect inputs');
end

syms z;

A = eye(2);
for i = 1:Nr,
   Ai = [1 rf(i);rf(i)*z z] * A;
   A  = Ai;
end

Hb = A(1,2);
Ha = A(2,2);
Gb = prod(tf)*(z^(Nr-2));
Ga = A(2,2);


return