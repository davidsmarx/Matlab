% get list of glasses from glass catalog
% glasscat = 'H:\ZEMAX lenses\Glasscat\OHARA.AGF';
glasscat = 'H:\ZEMAX lenses\Glasscat\SCHOTT_2000.AGF';

fid = fopen(glasscat,'rt');
cnt = 0;
while ~feof(fid),
   tline = fgetl(fid);
   if strcmp(tline(1:2),'NM'),
      cnt = cnt + 1;
      nsp = findstr(tline,' ');
      glass{cnt} = tline(nsp(1)+1:nsp(2)-1);
      formulacode{cnt} = sscanf(tline(nsp(2)+1:end),'%d',1);
      
      edline = fgetl(fid);
      if ~strcmp(edline(1:2),'ED'), error('ED line error'); end
      cte(cnt) = sscanf(edline(4:end),'%f',1);
      
      cdline = fgetl(fid);
      if ~strcmp(cdline(1:2),'CD'), error('CD line error'); end
      cd{cnt} = sscanf(cdline(3:end),'%f',[1 inf]);
      
      tdline = fgetl(fid);
      if ~strcmp(tdline(1:2),'TD'), error('TD line error'); end
      td{cnt} = sscanf(tdline(3:end),'%f',[1 inf]);
      
      odline = fgetl(fid);
      
      ldline = fgetl(fid);
      if ~strcmp(ldline(1:2),'LD'), error('LD line error'); end
      lamminmax{cnt} = sscanf(ldline(4:end),'%f %f');
   end
end
fclose(fid);

T = [0; 70];
for  ii = 1:cnt,
   switch formulacode{ii}
   case 1, formula = 'schott';
   case 2, formula = 'sellmeier1';
   otherwise error(['unknown formula code: ' num2str(formulacode{ii})]);
   end
   % index at room temp
   nrm(ii) = dispersionformula2index(cd{ii},formula,1.55,td{ii},20,23);
   
   % estimate dn/dT
   nn = dispersionformula2index(cd{ii},formula,1.55,td{ii},20,T);
   dn(ii) = 1e6*diff(nn)/diff(T);
end

figure, plot(cte,dn,'o'), grid, xlabel('CTE [x10^{-6}]'), ylabel('dn/dT')

strain = dn + nrm.*cte;
figure, plot(nrm,strain,'o'), grid, xlabel('index'), ylabel('thermo-optic strain [10^{-6}/\circC]')

% search for smallest thermo-optic strain
[dmin, imin] = sort(strain);
for ii = 1:10,
   ij = imin(ii);
   fprintf('%7s n = %5.2f cte = %5.2e dn = %6.3e strain = %e\n',...
      glass{ij},nrm(ij),cte(ij),dn(ij),strain(ij));
end

return

[reldispsort, isort] = sort(reldisp);
figure, plot(reldispsort), grid

isf57 = strmatch('SF57',glass,'exact');
disp([glass{isf57} ': relative dispersion = ' num2str(reldisp(isf57))]);

for n = cnt-20:cnt,
   disp([glass{isort(n)} ': rel disp = ' num2str(reldisp(isort(n)))...
         ' index = ' num2str(index(1,isort(n)))]);
end

return
