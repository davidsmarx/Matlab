function B = MyImopen(A, k)
% B = MyImopen(A, k)
% binary only at this time

Btmp = MyErosion(A, k);
B    = MyDilation(Btmp, k);

end % main


function B = MyErosion(A, k)

Btmp = conv2(double(A), k, 'same');

B = Btmp >= sum(k(:));

end % erosion

function B = MyDilation(A, k)

Btmp = conv2(double(A), k, 'same');

B = Btmp >= 1;

end % dilation
