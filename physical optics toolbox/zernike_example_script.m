x = linspace(-1,1)';
y = linspace(-1,1)';
[X, Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);

% non-zero Zernike 8 term:
Z = [zeros(1,7) 0.5];

% calculate Zernike function on non-rectangular grid:
phi = zeros(size(X));
phi(R<=1) = zernikeval(Z, X(R<=1), Y(R<=1), 1.0, 'Noll');

figure, imagesc(x,y,phi), axis image



