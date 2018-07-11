function z = besselj_zeros(nu,n)
% z = besselj_zeros(nu)
% z = besselj_zeros(nu,n)
%
% return the zeros of the besselJ_nu function corresponding to the integers
% in n. If no input arguments, then a vector of the first 1000 zeros is
% returned.
%
% nu = 0 to 10
%
% The tabulated data is stored in 'besselj[nu]_zeros.mat' and was produced by
% Mathematica with the following commands:
%
% In[1]:= << NumericalMath`BesselZeros`
% 
% In[2]:= z = BesselJZeros[0, 1000, AccuracyGoal -> 24, WorkingPrecision -> 34];
% 
% In[3]:= Export["C:\Documents and Settings\David\Matlab\physical optics \
% toolbox\besselj0_zeros.mat", z, "MAT"]

NUMAX = 10;

persistent zdata;

if isempty(zdata),
    for nu = 0 : NUMAX,
        s = load(['besselj' num2str(nu) '_zeros']);
        fn = fieldnames(s);
        ztmp = getfield(s,fn{1});
        zdata(:,nu+1) = ztmp(:);
        %zdata = zdata(:);
    end
end

switch nargin,
    case 1,
        if nu < 0 || nu > NUMAX || nu ~= round(nu),
            error(['nu = ' num2str(nu) ' is not supported']);
        end
        z = zdata(:,nu+1);

    case 2,
        if nu < 0 || nu > NUMAX || nu ~= round(nu),
            error(['nu = ' num2str(nu) ' is not supported']);
        end
        z = zdata(n,nu+1);
        
    otherwise,        
        error('usage: besselj_zeros(nu,n)');
end

return