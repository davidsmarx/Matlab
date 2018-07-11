n0 = 1;
n2 = 1.6;

t0 = linspace(0,89,90)'*P;

ns = (n0.*sin(t0)).^2;

n1s = sqrt( n0.*cos(t0).*sqrt(n2.^2 - ns) + ns );

t2 = asin((n0/n2)*sin(t0));
ap = n0*n2./(cos(t0).*cos(t2));
n1p = sqrt( 0.5*( ap + sqrt(ap.*ap - 4*ap.*ns)) );

figure, plot(t0/P,[n1s n1p]), legend('S','P'), grid


return

% one layer
syms z r1 r2; A = [z r2;r2*z 1]*[z r1;r1*z 1];
hb = subs(A(2,1),z,-1);
syms k0 k1 k2;
hn = subs(hb,{r1,r2},{(k1-k0)/(k1+k0),(k2-k1)/(k2+k1)});
hn = simplify(hn);
pretty(hn)

% two layers
syms z r1 r2 r3; A = [z r3;r3*z 1]*[z r2;r2*z 1]*[z r1;r1*z 1];
hb = subs(A(2,1),z,-1);
syms k0 k1 k2 k3;
hn = subs(hb,{r1,r2,r3},{(k1-k0)/(k1+k0),(k2-k1)/(k2+k1),(k3-k2)/(k3+k2)});
hn = simplify(hn);
pretty(hn)

return

n0 = 1;
n3 = 2;

t0 = linspace(0.01,89,90)'*P;
t3 = asin((n0/n3)*sin(t0));

ns = (n0.*sin(t0)).^2;

n1 = sqrt( ns.*(n0/n3 - 1)./(cos(t0)./cos(t3) - 1) );
n2 = n3*(n1.^2)/n0;

figure, plot(t0/P,[n1 n2]), legend('n1','n2'), grid
