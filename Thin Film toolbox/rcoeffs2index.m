function n = rcoeffs2index(r,nmed)
% n = rcoeffs2index(r,nmed)
% r = reflection coefficients at each boundary
% r(1) is at substrate, r(end) is at incident medium/first layer
% nmed = index of incident medium (default = 1)

% check reflection coefficients have magnitude < 1
if any(abs(r(2:end)) > 1.0), error('reflection coeffs with mag > 1 are not physical'); end

N = length(r);
n = zeros(1,N+1);

try n(end) = nmed; catch n(end) = 1; end

ii = r~=-1;
rr(ii) = (1 - r(ii))./(1 + r(ii));
rr(~ii) = inf;

for m = N:-1:1,
   n(m) = n(m+1)*rr(m);
end

return