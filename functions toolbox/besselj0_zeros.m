function z = besselj0_zeros(n)
% z = besselj0_zeros(n)
%
% return the zeros of the besselJ0 function corresponding to the integers
% in n. If no input arguments, then a vector of the first 1000 zeros is
% returned.
%
% The tabulated data is stored in 'besselj0_zeros.mat' and was produced by
% Mathematica with the following commands:
%
% In[1]:= << NumericalMath`BesselZeros`
% 
% In[2]:= z = BesselJZeros[0, 1000, AccuracyGoal -> 24, WorkingPrecision -> 34];
% 
% In[3]:= Export["C:\Documents and Settings\David\Matlab\physical optics \
% toolbox\besselj0_zeros.mat", z, "MAT"]


persistent zdata;

if isempty(zdata),
    s = load('besselj0_zeros');
    fn = fieldnames(s);
    zdata = getfield(s,fn{1});
    zdata = zdata(:);
end

if nargin == 0,
    z = zdata;
else,
    z = zdata(n);
end

return