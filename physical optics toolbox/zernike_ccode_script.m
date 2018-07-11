% sym r;
% maple('with','orthopoly')
% n = 3; m = 1;
% r.^m * maple('P',(n-m)/2,m+1,m+1,r^2)
% 
% Rnm = maple('sum(''(-1)^s * (n-s)! * r^(n-2*s) / (s! *((n+m)/2-s)! * ((n-m)/2-s)!)'',''s''=0..(n-m)/2)')
% 
% maple(['eval(' Rnm ', {n = ' num2str(n) ', m = ' num2str(m) '})']);

%%
syms n m r s t;
nsf = sym('(n-s)!');
sf = sym('s!');
npmf = sym('((n+m)/2-s)!');
nmmf = sym('((n-m)/2-s)!');
Rnm = maple('sum',(-1)^s * nsf * r^(n-2*s) / (sf * npmf * nmmf),'s=0..(n-m)/2');
cnt = 0;
for in = 0:10,
    for im = in:-2:0,
        ac = maple('eval',Rnm*cos(im*t),['{n = ' num2str(in) ', m = ' num2str(im) '}']);
        disp([num2str(cnt) ': ' num2str(in) ', ' num2str(im) ': ' ccode(maple('simplify',ac))]);
        cnt = cnt+1;
        
        as = maple('eval',Rnm*sin(im*t),['{n = ' num2str(in) ', m = ' num2str(im) '}']);
        if as ~= 0,
            disp([num2str(cnt) ': ' num2str(in) ', ' num2str(im) ': ' ccode(maple('simplify',as))]);
            cnt = cnt + 1;
        end
    end
end
