function [ii, distance] = findnearest(val, A)
% i = findnearest(val, A)
%
% like find, but for real numbers where find() doesn't work because real
% numbers are not exactly equal

[distance, ii] = min(abs(A - val));
