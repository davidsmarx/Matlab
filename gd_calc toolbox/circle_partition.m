function x=circle_partition(N)
% x=circle_partition(N)
%
% Partition the first quadrant of the unit circle into N rectangular blocks
% with corner vertices at [x(1),x(end)], [x(2),x(end-1)], ...
% [x(end),x(1)]. x(j) is monotonic decreasing with j.
%
% Documentation reference: GD-Calc_Demo.pdf, Appendix B.
%
% Version 11/05/2005
% Author: Kenneth C. Johnson
% software.kjinnovation.com
% This file is Public Domain software.

theta1L=(pi/4)/(sqrt(3)*(N-1/2)+1)^(2/3);
theta1H=(pi/4)/(sqrt(3)*(N+1/2)+1)^(2/3);
errL=target_err(theta1L,N);
errH=target_err(theta1H,N);
theta1=theta1L+(theta1H-theta1L)*(-errL)/(errH-errL);
d_theta1=theta1H-theta1L;
while abs(errH-errL)>sqrt(eps)
    if errL*errH<=0
        d_theta1=d_theta1/10;
    end
    theta1L=theta1-d_theta1/2;
    theta1H=theta1+d_theta1/2;
    errL=target_err(theta1L,N);
    errH=target_err(theta1H,N);
    theta1=theta1L+(theta1H-theta1L)*(-errL)/(errH-errL);
end
area=theta1-0.5*sin(2*theta1);
theta=[theta1,zeros(1,N-1)];
outside_pt=true;
for j=2:N
    theta(j)=next(theta(j-1),area,outside_pt);
    outside_pt=~outside_pt;
end
theta_=next(theta(end),area,outside_pt);
theta=[theta,pi/2-theta(end:-1:1)];
x=cos(theta(1:2:end));

function theta=next(theta,area,outside_pt)
theta_=theta+2*sqrt(area/sin(2*theta));
prev_diff=[];
c=cos(theta);
s=sin(theta);
while true
    c_=cos(theta_);
    s_=sin(theta_);
    if outside_pt
        theta_=theta_-(c*s_-0.5*(c*s+c_*s_)-0.5*(theta_-theta)-area)/...
            (c*c_-0.5*(c_*c_-s_*s_)-0.5);
    else
        theta_=theta_-(c_*s-0.5*(c*s+c_*s_)+0.5*(theta_-theta)-area)/...
            (-s_*s-0.5*(c_*c_-s_*s_)+0.5);
    end
    if isempty(prev_diff)
        prev_diff=theta_-theta;
    else
        diff=theta_-prev_theta_;
        if abs(diff)<eps || abs(diff)>=abs(prev_diff)
            break
        end
        prev_diff=diff;
    end
    prev_theta_=theta_;
end
theta=theta_;

function err=target_err(theta1,N)
area=theta1-0.5*sin(2*theta1);
theta=theta1;
outside_pt=true;
for j=2:N
    theta=next(theta,area,outside_pt);
    outside_pt=~outside_pt;
end
theta_=next(theta,area,outside_pt);
err=theta+theta_-pi/2;
