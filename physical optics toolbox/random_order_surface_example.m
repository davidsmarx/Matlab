% xx = random order list of x-values
% yy = random order list of y-values
%  z = list of z-values for each xx,yy
A = textread('data_example.txt');
xx = A(:,1)/14;
yy = A(:,2)/14;
 z = A(:,3);

stepsize = 1;
x = [min(xx):stepsize:max(xx)];
y = [min(yy):stepsize:max(yy)];
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));

for ii = 1:length(X(:)),
    atmp = z(xx==X(ii) & yy==Y(ii));
    if ~isempty(atmp)
	    Z(ii) = atmp;
    end
end

figure, pcolor(X,Y,Z)
figure, imagesc(x,y,Z)

% pcolor default:
set(gca,'YDir','normal')

% imagesc default:
set(gca,'YDir','reverse')