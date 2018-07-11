function [Surfs, zchan] = z_Getlensdata(zchan,nconfig)
% [Surfs, zchan] = z_Getlensdata(zchan,nconfig)
% alternative to z_Getlens()
% nconfig (optional, default = 1), get lens data for configuratin #nconfig
% zchan (optional, default = initiate new lens) previously initiated channel
% Surfs is array of surface structs
% each surface is a struct with the following fields:
% nsurf, type, curv, thick, glass, radius, conic, parms, coat, name

if nargout == 0 & nargin == 0,
   error('usage: [Surfs, zchan] = z_Getlensdata(nconfig,zchan)');
end

% check if zchan exists and is valid, and get system data
try,
   sysstr = ddereq(zchan,'GetSystem',[1 1]);
catch,
   % initiate contact with ZEMAX and load lens data into ZEMAX buffer
   zchan = ddeinit('zemax','topic');
   retstr = ddereq(zchan,'GetRefresh',[1 1]);
   sysstr = ddereq(zchan,'GetSystem',[1 1]);
end

% get # of surfaces
nsurf = sscanf(sysstr,'%d',1);

% set the configuration
try,
   errcond = z_setconfig(zchan,nconfig);
end
[cf, ncf] = z_getconfig(zchan);
disp(['current configuration: #' num2str(cf) ' out of ' num2str(ncf)]);

% get each surface and store as a struct (surface 0 is object)
for n = 1:nsurf,
   % get the data for each surface and store as struct
   surfst = z_getsurf(zchan,n);

   % get names of the surface
   surflab = deblank( ddereq(zchan,['GetComment,' num2str(n)],[1 1]) );
   if isempty(surflab),
      surflab = ['surf_' num2str(n)];
   end
   % replace spaces with _
   ii = findstr(surflab,' '); surflab(ii) = '_';
   
   surfst.name = surflab;
   Surfs(n) = surfst;
   
end

return


