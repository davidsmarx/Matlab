% script to apply scallop correction and then focus with seismic
% the script asks for the following user inputs:
%   base name of view-formatted sonar data
%   name of the 'imswat.in' file used to produce the sonar data
%   name for the view-formatted result

speedc = 1500; % [m/s]

basename = input('base name of view-formatted data: ');
[cdata parms] = sasview(basename);
Ny = parms(5); ymin = parms(6); dy = parms(7);

imswatdata = input('name of imswat .in file: ');
if imswatdata(length(imswatdata)-2:length(imswatdata)) ~= '.in',
   imswatdata = [imswatdata '.in'];
end
if exist(imswatdata) ~= 2, error([imswatdata 'does not exist']); end

tic;

imswat = load(imswatdata);
                              % imswat(15) is distance to target [m]
                              % imswat(16) is image width [m] (range direction)
                              % imswat(17) is image height [m] (az direction)
freq     = imswat(23)*1000;   % imswat(23) is the frequency in kHz
N_el     = imswat(27);        % imswat(27) is the number of receiver elements
dx_l     = imswat(28);        % imswat(28) is the element separation
off_tr   = imswat(53);        % imswat(53) is the projector offset
velocity = imswat(54)*0.5144; % imswat(54) is the platform velocity in knots
fprintf('freq = %f Hz\n',freq);
fprintf('array: %d elements spaced %fm\n',N_el,dx_l);
fprintf('projector offset: %f m\n',off_tr);
fprintf('platform velocity = %f m/s\n',velocity);

y = 1; n = 0;
doscallop = input('does this data require scallop correction? ');
if doscallop,
   t = (ymin+dy*[0:Ny-1])*2/speedc;
   x_el = ([-0.5*(N_el-1):0.5*(N_el-1)])*dx_l;
   
   rcs = scallop_rshift(cdata,t,x_el,off_tr,velocity,freq);
%   rcs = scallop(cdata,t,x_el,x_tran,velocity,freq);

   disp('scallop correction complete, now doing seismic...');
else,
   rcs = cdata;
end
clear cdata;  % got to save precious memory

outname = [basename '_sas'];
tmp = input(['base name of view-formatted focused result [' outname ']: ']);
if ~isempty(tmp), outname = tmp; end

fid = fopen(outname,'w');
while fid == -1,
  outname = input(['cannot open file ' outname ' new basename: ']);
  fid = fopen(outname,'w');
end
fclose(fid);

if computer == 'PCWIN'
   crec = seismic(rcs,parms,parms(6),speedc/freq,0,0);
   cnt = saswrite(outname,crec,parms);
   
else,
   scallopname = ['/home/dsm/cssmocomp/tmpscallop'];
   cnt = saswrite(scallopname,rcs,parms);

   disp(['! nice seismic ' scallopname ' ' outname ' ' num2str(parms(6)) ' ' num2str(1500/freq) ' 0 0']);
   eval(['! nice seismic ' scallopname ' ' outname ' ' num2str(parms(6)) ' ' num2str(1500/freq) ' 0 0']);
   [crec, prec] = sasview(outname);
end

disp('variable ''crec'' is the result focused image');
fprintf('processing time: %.2f minutes\n',toc/60);

return
