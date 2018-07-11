function gd = read_schott_cat(gname)
% gd = read_schott_cat(gname)
%
% gname (optional) is a string glass name, 'all' to get the whole
% catalog, help to launch the online help for catalog explanations
% found in "C:\Program Files\SCHOTT\OGLCatalogue2000\SCHOTT2000E.HLP"
%
% gd is a struct contianing all the catalog parameters
% colorCode is not read correctly
% 

if nargin == 0, gname = 'all'; end

if strcmpi(gname,'help'),
    dos('C:\Program Files\SCHOTT\OGLCatalogue2000\SCHOTT2000E.HLP &');
    return
end

fid = fopen('Schott Data 2000.asc','rt');

% if gname is 'all', read all the glass types in an array of structs
if strcmpi(gname,'all'),
   cnt = 0;
   tline = fgetl(fid);
   while tline ~= -1,
      cnt = cnt+1;
      gd(cnt) = read_data(fid,tline);
%       gd(cnt) = read_data(fid,gd(cnt));
      tline = fgetl(fid);
   end
   
% else search for the desired glass type
else,
   gd = read_data(fid,goto_glass(fid,gname));
   
end

fclose(fid);

return

function nameout = goto_glass(fid,gname)

while 1,
   tline = fgetl(fid);
   if strcmpi(tline,gname), break, end
   if ~isempty(tline) & tline == -1, error(['glass ' gname ' not found']); end
end

nameout = tline;

return


function gd = read_data(fid,gname)
% read the parameters up to the next blank line

gd.name = gname;

tline = fgetl(fid); 
while ~isempty(tline),
   try,
      [fld, val] = strread(tline,'%s%f','delimiter',':');
      fldstr = fld{1}; 
      fldstr(findstr(fldstr,'''')) = 'p';
      fldstr(findstr(fldstr,'-')) = 'm';
   catch,
      [fld, valc] = strread(tline,'%s%s','delimiter',':');
      fldstr = fld{1};
      val = valc{1};
   end
   
   if isempty(val), val = 0; end
   gd = setfield(gd,fldstr,val);
   tline = fgetl(fid); 
end

% convert constants of dispersion formula to single polynomial
gd.disp_poly = [gd.B1 gd.C1 gd.B2 gd.C2 gd.B3 gd.C3];
% convert constants of dn/dT to single polynomial
gd.dndt_poly = [gd.D0 gd.D1 gd.D2 gd.E0 gd.E1 gd.lamdaTK];

return