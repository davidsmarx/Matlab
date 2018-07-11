function yi = interpcircle(theta, y, ti, method)
% yi = interpcircle(theta, y, ti)
% yi = interpcircle(theta, y, ti, 'method')
%
% interpcircle calls interp1 with method 'method'
%
% assumes theta is monotonic ascending or descending and
% |theta(end) - theta(1)| <= 2 pi

if ~exist('method','var') | isempty(method),
    method = [];
end

% check if end points need to be wrapped
if abs(theta(end) - theta(1)) < 2*pi,

    %theta = unwrap(theta);
    
    % is theta ascending or descending
    if theta(end) > theta(1),
        % ascending
        tt = [theta(end) - 2*pi; theta(:); theta(1) + 2*pi];
        yy = [y(end); y(:); y(1)];
    else
        % descending
        tt = [theta(end) + 2*pi; theta; theta(1) - 2*pi];
        yy = [y(end); y(:); y(1)];
    end

else
    % no need to wrap, interp as is
    tt = theta;
    yy = y;
    
end

yi = interp1(tt, yy, ti, method);

