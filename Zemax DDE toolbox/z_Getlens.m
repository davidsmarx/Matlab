function [Lens, zchan] = z_Getlens(zchan,nconfig)
% [Lens, zchan] = z_Getlens
% Lens = z_Getlens(zchan, nconfig)
% [Lens, zchan] = z_Getlens([], nconfig)
%
% nconfig (optional, default = 1), get lens data for configuratin #nconfig
% zchan (optional, default = initiate new lens) previously initiated channel
% Lens is struct with a field for each surface
% each surface is a struct with the following fields:
% nsurf, type, curv, thick, glass, radius, conic, parms, coat

% parse inputs
if ~exist('nconfig','var') | isempty(nconfig), nconfig = 1; end
if exist('zchan','var') & ~isempty(zchan),
    % check if zchan exists and is valid, and get system data
    try,
        sysstr = ddereq(zchan,'GetSystem',[1 1]);
        validzchan = true;
    catch,
        % error occured
        validzchan = false;
    end
else
    validzchan = false;
end

% if necessary, initiate contact with ZEMAX and load lens data into ZEMAX
% buffer
if validzchan == false,
    zchan = ddeinit('zemax','topic');
    retstr = ddereq(zchan,'GetRefresh',[1 1]);
    sysstr = ddereq(zchan,'GetSystem',[1 1]);
end

% get # of surfaces
nsurf = sscanf(sysstr,'%d',1);

% set the configuration
z_setconfig(zchan,nconfig);
[cf, ncf] = z_getconfig(zchan);
%disp(['current configuration: #' num2str(cf) ' out of ' num2str(ncf)]);

Lens = struct;
% get each surface and store as a struct
for n = 0:nsurf,
   % get the data for each surface and store as struct
   surfst = z_getsurf(zchan,n);

   % get names of the surfaces and use as variable names
   surflab = deblank( ddereq(zchan,['GetComment,' num2str(n)],[1 1]) );
   if isempty(surflab),
      surflab = ['surf_' num2str(n)];
   end
   % replace invalid characters with _
   surflab(regexp(surflab,'\W')) = '_';
   % cannot start with a number
   if regexp(surflab(1),'\d'), surflab = ['s_' surflab]; end
   
   % check for identical surface labels
   if isfield(Lens,surflab),
      surflab = [surflab '_' num2str(n)];
      %disp(['Warning: identical surface name, changing to: ' surflab]);
   end
   
   % add to the Lens struct
   try,
       Lens.(surflab) = surfst;
       %eval(['Lens.' surflab '= surfst;']);
   catch
      disp('z_Getlens: error in eval');
      keyboard;
   end
   
end

%disp(fieldnames(Lens));

return


