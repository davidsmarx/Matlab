function [gd, glasscat] = read_zemax_glasscat(glasscatin,gname)
% [gd, glasscat] = read_zemax_glasscat
% gd = read_zemax_glasscat(glasscat)
% gd = read_zemax_glasscat(glasscat,gname)
%
% example:
%    gd = read_zemax_glasscat('INFRARED.AGF','SILICON');

global UM MM;

if isdir('C:\Program Files\ZEMAX\Glasscat\'),
    ZEMAXGLASSCATPATH = 'C:\Program Files\ZEMAX\Glasscat\';
elseif isdir('C:\ZEMAX\Glasscat\'),
    ZEMAXGLASSCATPATH = 'C:\ZEMAX\Glasscat\';
elseif isdir('C:\Documents and Settings\All Users\Application Data\ZEMAX\Glasscat\'),
    ZEMAXGLASSCATPATH = 'C:\Documents and Settings\All Users\Application Data\ZEMAX\Glasscat\';
else
    error('where is the ZEMAX glass catalog?');
end

if ~exist('glasscatin','var') | isempty(glasscatin),
    [glasscat, pathname, fi] =...
        uigetfile({'*.agf'},'Select Glass Catalog',ZEMAXGLASSCATPATH);
    if isequal(glasscat,0), return, end
    glasscat = [pathname glasscat];
else,
    glasscat = [ZEMAXGLASSCATPATH glasscatin];
end

fid = fopen(glasscat,'rt');
if fid == -1, error(['error opening ' glasscat]); end

%try,
cnt = 0;
while ~feof(fid),
   tline = fgetl(fid);
   switch tline(1:2),
       
       case 'NM',
           cnt = cnt + 1;
           nsp = findstr(tline,' ');
           gd(cnt).name = tline(nsp(1)+1:nsp(2)-1);
           ctmp = textscan(tline(nsp(2)+1:end),'%f','TreatAsEmpty','-');
           atmp = ctmp{1};
           
           gd(cnt).formulacode = atmp(1);
           switch gd(cnt).formulacode,
               case 1, gd(cnt).formula = 'schott';
               case 2, gd(cnt).formula = 'sellmeier1';
               case 3, gd(cnt).formula = 'herzberger';
               case 5, gd(cnt).formula = 'conrady';
               case 6, gd(cnt).formula = 'sellmeier3';
               case 7, gd(cnt).formula = 'handbook1';
               case 8, gd(cnt).formula = 'handbook2';
               case 9, gd(cnt).formula = 'sellmeier4';
               case 11, gd(cnt).formula = 'sellmeier5';
               otherwise warning([gd(cnt).name ': unknown formula code: ' num2str(gd(cnt).formulacode)]);
           end

           gd(cnt).nd = atmp(3);
           gd(cnt).vd = atmp(4);
           

           statuscode = atmp(6);
           switch statuscode,
               case 0, gd(cnt).status = 'standard';
               case 1, gd(cnt).status = 'preferred';
               case 2, gd(cnt).status = 'obsolete';
               case 3, gd(cnt).status = 'special';
               case 4, gd(cnt).status = 'melt';
               otherwise, gd(cnt).status = 'unknown';
           end
           
           if length(atmp) >= 7,
               gd(cnt).meltfreq = atmp(7);
           else
               gd(cnt).meltfreq = [];
           end
      
       case 'GC',
           % don't know what this is
           
       case 'ED',
           A = textscan(tline,'ED %f %f %f %f %f','treatAsEmpty','-');
           gd(cnt).cte = A{1};
           gd(cnt).p = A{3};
           gd(cnt).dPgF = A{4};
      
       case 'CD',
           gd(cnt).disp_poly = sscanf(tline(3:end),'%f',[1 inf]);
      
       case 'TD',
           gd(cnt).dndt_poly = sscanf(tline(3:end),'%f',[1 inf]);
      
       case 'OD',
           A = textscan(tline,'OD %f %f %f %f %f %f','treatAsEmpty','-');
           if tline(4) == '-',
               atmp = sscanf(tline(5:end),'%f',[1 inf]);
               atmp = [-1 atmp];
           else,
               atmp = sscanf(tline(4:end),'%f',[1 inf]);
           end

           gd(cnt).RelCost = A{1};
           gd(cnt).CR = A{2};
           gd(cnt).FR = A{3};
           gd(cnt).SR = A{4};
           gd(cnt).AR = A{5};
           gd(cnt).PR = A{6};
      
       case 'LD',
           atmp = sscanf(tline(4:end),'%f %f');
           gd(cnt).MinWave = atmp(1)*UM;
           gd(cnt).MaxWave = atmp(2)*UM;
           
       case 'IT',
           % not yet implemented
           
       otherwise,
           % don't know this code
   end

end

% catch,
%     fclose(fid);
%     error(lasterr);
% end

fclose(fid);

if exist('gname','var') & ~isempty(gname),
    ii = strmatch(gname,{gd.name},'exact');
    if isempty(ii),
        error(['glass ' gname ' not found']);
    end
    gd = gd(ii);
end

