lam = [0.95 1.375 1.7]*UM; % (micron)

%glass = read_schott_cat;
[glass, glasscat] = read_zemax_glasscat;
%glasscatlist = dir('C:\zemax\glasscat\*.agf');

% for ii = 1:length(glasscatlist),
%     glasscat = glasscatlist(ii).name;
%     glass = read_zemax_glasscat(glasscat);
    
v = zeros(size(glass));
n1300 = zeros(size(glass));
for ii = 1:length(v),
    
    if lam(1) > glass(ii).MinWave & lam(end) < glass(ii).MaxWave ...
            & glass(ii).RelCost > 0,
        n = dispersionformula2index(glass(ii).disp_poly,glass(ii).formula,lam);
        if ~any(isnan(n)),
            n1300(ii) = n(2);
            v(ii) = (n(2) - 1)./(n(1) - n(3));
        else
            n1300(ii)= NaN;
            v(ii) = NaN;
        end
    else
        n1300(ii) = NaN;
        v(ii) = NaN;
    end
   
end

figure, plot(n1300,v,'.'), grid,...
   text(n1300,v,{glass.name},'fontsize',12),...
   xlabel('index of refraction, n1300'), ylabel('Abbe Number v'),...
   title(pwd2titlestr(glasscat))

% pause;
% 
% end
