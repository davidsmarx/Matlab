function rmsresidual = cumulativeOPDresidual(x,y,opd,nmax,R)
% rmsresidual = cumulativeOPDresidual(x,y,opd,nmax,R)

global PM;

rmsresidual = zeros(nmax,1);

for ii = 1:nmax,
    [Z, res, rmsresidual(ii)] = zernikefit(x,y,opd,ii,R);
end

if nargout == 0,
    figure, plot([2:nmax],rmsresidual(2:end)/PM,'-o'), grid
end
