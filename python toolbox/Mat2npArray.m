function nparray = Mat2npArray(A)
% nparray = Mat2npArray(A)
%
% a scalar is translated to a numpy scalar (shape = ())
% if you want to translate a scalar to numpy array shape = (1,)
% do this:
%    py.numpy.reshape(py.numpy.array(a), 1);

if isscalar(A),
    nparray = A;
    return;
end

% rearrange dimensions
szA = size(A);
nd = ndims(A);
if nd > 2
    At = permute(A, [2, 1, 3:nd]);
    newsize = [szA(3:nd) szA(1) szA(2)];
elseif nd == 2,
    At = A.';
    newsize = size(A);
    % if a row vector, make a numpy array of size N,
    if newsize(1) == 1, newsize = newsize(2); end
else,
    % should never get here
    error('unknown A dims');
        
end

% template for matlab array to numpy array
fMat2Py = @(At, pytype, newsize) py.numpy.reshape(py.numpy.array(At(:).', pytype), int32(newsize));


% check type
mattype = class(A);
switch mattype,
    case 'double',
        if isreal(A),
            nparray = fMat2Py(At, 'float64', newsize);

        else
            pyI = py.numpy.complex128(1i);
            pyAr = fMat2Py(real(At), 'complex128', newsize);
            pyAi = fMat2Py(imag(At), 'complex128', newsize);
            nparray = pyAr + pyI*pyAi;
        end
        
    case {'int32', 'int8', 'int16', 'int64'},
        nparray = fMat2Py(At, mattype, newsize);
        
    case 'logical',
        At = zeros(size(A));
        At(A) = 1;
        nparray = fMat2Py(At, 'bool', newsize);

    otherwise
        error(['mat type ' mattype ' not handled yet']);
end

% re-transpose back
if nd > 2,
    nparray = py.numpy.transpose(nparray, int32([nd-1, nd, 1:nd-2]-1));
end

