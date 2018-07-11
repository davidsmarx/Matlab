unitsdefinitions;

h = actxserver('int_met_comserver.Diffraction')

[ A, dx, dy, lambda, curv, x, y ] = ReadWavefrontUNF( 'ns_applyfieldsep.unf' );
figure, imagesc(x,y,abs(A)), axis image

set(h,'wavelength',lambda)
set(h,'dx',dx);
set(h,'dy',dy);
set(h,'curvature',curv);
[c, d] = invoke(h,'focuslens',real(A),imag(A));
B = c+j*d;
lambda = get(h,'wavelength');
dx = get(h,'dx');
dy = get(h,'dy');
curv = get(h,'curv');
xx = invoke(h,'x_vector_get',x);
yy = invoke(h,'y_vector_get',y);
figure, imagesc(xx,yy,abs(B)), axis image

delete(h);