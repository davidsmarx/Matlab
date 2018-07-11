function y = myfftshift(x)
% y = myfftshift(x)
% if x is a vector, y = fftshift(x)
% if x is a matrix, columns of y = fftshift(columns(x))
% Note: myfftshift is different from fftshift for matrices
% because fftshift(x=matrix) interchanges quadrants

[nr nc] = size(x);

if ( (nr==1) | (nc==1) ),
   y = fftshift(x);
   return;
end

nr2 = floor(nr/2);

y = zeros([nr nc]);

%fprintf('size of y = %d x %d\n',size(y));

y(1:nr2,:) = x(nr-nr2+1:nr,:);
y(nr2+1:nr,:) = x(1:(nr-nr2),:);
