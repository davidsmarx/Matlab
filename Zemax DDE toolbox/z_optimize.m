function val = z_optimize(zchan,n,timeout)
% val = z_optimize(zchan,n,timeout)
% n: number of cycles
% where n is the number of cycles to run. The return value is the final merit function. If the merit function value
% returned is 9.0E+009, the optimization failed, usually because the lens or merit function could not be evaluated.
% If n is zero, the optimization runs in automatic mode. If n is less than zero (for example, n = -1), Optimize updates
% all operands in the merit function and returns the current merit
% function, and no optimization is performed.
% 
% timeout (optional) = timeout for waiting for response (s), default = 30s

global MS;

if ~exist('timeout','var')
    timeout = 30;
end

cmdstr = ['Optimize,' num2str(n)];

retstr = ddereq(zchan,cmdstr,[1 1],timeout/MS);
if retstr == 0, warning('timeout condition during optimization'); end
if ~ischar(retstr), warning('return format error'); keyboard; end

val = sscanf(retstr,'%f,',[1 inf]);

return


