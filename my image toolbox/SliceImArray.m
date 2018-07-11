function imout = SliceImArray(imin, N)
% imout = SliceImArray(imin, [Ny Nx])

if any(N > size(imin)),
    imout = imin;
    return
end

[xi, yi] = CreateGrid(imin);

imout = zeros(N);
[xo, yo] = CreateGrid(imout);


imout = imin(yi >= yo(1) & yi <= yo(end), xi >= xo(1) & xi <= xo(end));
