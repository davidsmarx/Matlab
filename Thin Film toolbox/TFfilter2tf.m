function [HB, HA, GB, GA, rf, rb, tf, tb] = TFfilter2tf(n)
% [Hb, Ha, Gb, Ga, r] = TFfilter2tf(n)
% n(1) = index of substrate
% n(2:end-1) = indexes of QWOT
% n(end) = index of incident medium
% Hb = numerator polynomial of z-transform of reflected field
% Ha = denominator polynomial of z-transform of reflected field
% Gb = numerator polynomial of z-transform of transmitted field
% Ga = denominator polynomial of z-transform of transmitted field
% rf = reflection coefficients at each thin film boundary
% rb = -rf (reflection coefficient from other direction)
% tf = forward transmission coeff. from n(2) to n(1) i.e. incident med. to sub.
% tb = backward transmission coeff. from n(1) to n(2) i.e. sub. to inc. med.
%
% This function first calculates the reflection and transmission
% coefficients r and t at each boundary, then calls rcoeffs2tf

% first calculate reflection coefficients
n0 = n(1:end-1);
n1 = n(2:end);
rf = (n1 - n0)./(n1 + n0);
rb = -rf;
tf = 2*n0./(n1 + n0);
tb = 2*n1./(n1 + n0);

% call rcoeffs2tf
%[HB, HA, GB, GA] = rcoeffs2tf(rf,rb,tf,tb);
[HB, HA, GB, GA] = rcoeffs2tf(rf);

return