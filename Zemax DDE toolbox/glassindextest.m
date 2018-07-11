% get list of glasses from glass catalog
glasscat = 'H:\ZEMAX lenses\Glasscat\SCHOTT.AGF';
fid = fopen(glasscat,'rt');
cnt = 0;
while ~feof(fid),
   tline = fgetl(fid);
   if strcmp(tline(1:2),'NM'),
      cnt = cnt + 1;
      nsp = findstr(tline,' ');
      glass{cnt} = tline(nsp(1)+1:nsp(2)-1);
      while ~strcmp(tline(1:2),'LD'),
         tline = fgetl(fid);
      end
      lamminmax{cnt} = sscanf(tline(4:end),'%f %f');
   end
end

% load lens
[Lens, zchan] = z_Getlens;
dumsf = Lens.LENS_FACET;

% load each glass into zemax and check indexes at defined wavelengths
for n = 1:cnt,
   if lamminmax{n}(1) < 1.5 & lamminmax{n}(2) > 1.61,
      dumsf.glass = glass{n};
      errcond = z_setsurf(zchan,dumsf);
      index(:,n) = z_getindex(zchan,dumsf);
      ndif(n) = index(2,n) - index(3,n);
      reldisp(n) = (index(1,n)-1)./ndif(n);
   end
end

[reldispsort, isort] = sort(reldisp);
figure, plot(reldispsort), grid

isf57 = strmatch('SF57',glass,'exact');
disp([glass{isf57} ': relative dispersion = ' num2str(reldisp(isf57))]);

for n = cnt-20:cnt,
   disp([glass{isort(n)} ': rel disp = ' num2str(reldisp(isort(n)))...
         ' index = ' num2str(index(1,isort(n)))]);
end

return
         
rc = ddeterm(zchan);