function [xf, yf, mag, theta] = ReadCodeVFMAdistortion(fn)
% [xf, yf, mag, theta] = ReadCodeVFMAdistortion(fn)
%
% xf, yf are same units as field definition in the Code V
% mag is same units as image plane in the Code V
% theta is (radians)

U = CConstants;

fid = fopen(fn,'rt');
if isequal(fid,-1), error(['cannot open file ' fn]); end

ltmp = fgetl(fid);
while isempty(strfind(ltmp, 'X-Field'))
    ltmp = fgetl(fid);
end
% data follows
A = textscan(fid,'%f %f %f %f', 101*101, 'delimiter', '\t');

fclose(fid);

[xf, yf, mag, theta] = deal(A{:});
theta = theta*U.P; % degrees to radians


