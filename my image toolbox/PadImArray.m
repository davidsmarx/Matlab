function [Imout, npad_pre, npad_post] = PadImArray(Im, N)
% Imout = PadImArray(Im, N)
%
% N is the size of the new array
% assumes N = [Ny Nx] >= size(Im)
%
% see also CropImage()

    Nim = size(Im);
    if any(N - Nim <= 0),
        %         warning('no padding necessary, returning Im unchanged');
        Imout = Im;
        return
    end

    
    npad_pre = ceil((N-Nim)/2);
    npad_post= floor((N-Nim)/2);
    
    % pdarray is in the image processing toolbox, and too often I can't get
    % a license from the license server for the image processing toolbox
    %     if all(npad_pre == npad_post),
    %         Imout = padarray(Im,npad_pre,'both');
    %     else
    %         Imout = padarray(Im,npad_pre,'pre');
    %         Imout = padarray(Imout,npad_post,'post');
    %     end
    
    Imtmp = padarray_pre(Im, npad_pre);
    Imout = padarray_post(Imtmp, npad_post);

end % main

function Aout = padarray_post(Ain, npad)

    Nin  = size(Ain);
    Nout = Nin + npad;

    Aout = cast(zeros(Nout),'like',Ain);
    
    Aout(1:Nin(1),1:Nin(2)) = Ain;

end

function Aout = padarray_pre(Ain, npad)

    Nin  = size(Ain);
    Nout = Nin + npad;
    
    Aout = cast(zeros(Nout),'like',Ain);

    Aout((npad(1)+1):end,(npad(2)+1):end) = Ain;

end
