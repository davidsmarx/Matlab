function B = dct2(A,isign)
%B = dct2(A,isign)
%2-dimensional discrete cosine transform
%A and B must be 2-d matrices with Nx and Ny powers of 2
%isign = 1 for forward transform and -1 for inverse transform

[nx ny] = size(A);
B = zeros([nx ny]);
Ap= zeros([nx ny]);

Ap = dct_(A,1);
B = dct_(Ap',1)';





if 0,
for i=1:ny,
   Ap(:,i) = dct(A(:,i));
end
for i=1:nx,
   B(i,:) = dct(Ap(i,:));
end
end