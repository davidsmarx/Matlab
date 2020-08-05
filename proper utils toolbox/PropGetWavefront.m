function [E, x, y, X, Y, R, T] = PropGetWavefront(wavefront, varargin)
% [E, x, y, X, Y, R] = PropGetWavefront(wavefront, varargin)
%
% combines Proper 'get' functions to return Efield and CreateGrid()

E = prop_get_wavefront(wavefront);
dx = prop_get_sampling(wavefront);

[x, y, X, Y, R, T] = CreateGrid(E, dx);
