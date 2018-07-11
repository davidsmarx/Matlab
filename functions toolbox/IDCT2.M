function B = idct2(A)
%B = dct2(A)
%2-dimensional inverse discrete cosine transform
%A and B must be 2-d matrices with Nx and Ny powers of 2

[nx ny] = size(A);
B = zeros([nx ny]);
Ap= zeros([nx ny]);

if 1,
for i=1:ny,
   Ap(:,i) = dct_(A(:,i),-1);
end
for i=1:nx,
   B(i,:) = dct_(Ap(i,:)',-1)';
end
end

if 0,
for i=1:ny,
   Ap(:,i) = idct(A(:,i));
end
for i=1:nx,
   B(i,:) = idct(Ap(i,:));
end
end
