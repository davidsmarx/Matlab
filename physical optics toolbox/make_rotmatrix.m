function R = make_rotmatrix(t_x,t_y,t_z)
% R = make_rotmatrix(t_x,t_y,t_z)
% R = make_rotmatrix([t_x t_y t_z])
%
% t_x, t_y, t_z are rotations about x, y, z-axis
% the order is R = Rz*Ry*Rx

switch nargin
    case 1,
        tmp = t_x;
        t_x = tmp(1);
        t_y = tmp(2);
        t_z = tmp(3);
    case 3,
        % do nothing
    otherwise,
        error('usage: R = make_rotmatrix(t_x,t_y,t_z)');
end

R = [
    cos(t_z) -sin(t_z)  0
    sin(t_z)  cos(t_z)  0
    0         0         1
    ] * [
    cos(t_y)  0         sin(t_y)
    0         1         0
   -sin(t_y)  0         cos(t_y)
    ] * [
    1         0         0
    0         cos(t_x) -sin(t_x)
    0         sin(t_x)  cos(t_x)
    ];


    