function result = NParray2Mat(nparray)
% A = NParray2Mat(nparray)
%
% taken from:
% https://www.mathworks.com/matlabcentral/answers/157347-convert-python-numpy-array-to-double

switch nparray.dtype.name.string
    
    case 'float64',       
        
        result = ConvertShape(nparray, @double, 'd');
        
    case 'complex128',
        rtmp = ConvertShape(nparray.real, @double, 'd');
        itmp = ConvertShape(nparray.imag, @double, 'd');
        result = rtmp + 1i*itmp;
        
    case {'int64','int32'}
        result = ConvertShape(nparray, @int64, 'l');
        
    otherwise,

        
end % switch type

end % main

function result = ConvertShape(nparray, funType, typechar)


        %nparray2mat Convert an nparray from numpy to a Matlab array
        %   Convert an n-dimensional nparray into an equivalent Matlab array
        data_size = cellfun(@int64,cell(nparray.shape));
        if length(data_size)==1
            % This is a simple operation
            result=funType(py.array.array(typechar, py.numpy.nditer(nparray)));
        elseif length(data_size)==2
            % order='F' is used to get data in column-major order (as in Fortran
            % 'F' and Matlab)
            result=reshape(funType(py.array.array(typechar, ...
                py.numpy.nditer(nparray, pyargs('order', 'F')))), ...
                data_size);
        else
            % For multidimensional arrays more manipulation is required
            % First recover in python order (C contiguous order)
            result=funType(py.array.array(typechar, ...
                py.numpy.nditer(nparray, pyargs('order', 'C'))));
            % Switch the order of the dimensions (as Python views this in the
            % opposite order to Matlab) and reshape to the corresponding C-like
            % array
            result=reshape(result,fliplr(data_size));
            % Now transpose rows and columns of the 2D sub-arrays to arrive at the
            % correct Matlab structuring
            result=permute(result,[length(data_size):-1:1]);
            
        end % if length(data_size)
        
        
end