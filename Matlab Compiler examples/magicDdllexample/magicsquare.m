function m = magicsquare(n)
%MAGICSQUARE generates a magic square matrix of the size specified
%    by the input parameter n.

% Copyright 2003-2005 The MathWorks, Inc.

if ischar(n)
    n=str2num(n);
end
m = magic(n)