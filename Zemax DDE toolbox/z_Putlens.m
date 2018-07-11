function zchan = z_Putlens(Lens,zchan)
% zchan = z_Putlens(Lens,zchan)
% Lens is a struct of surfaces as defined in z_Getlens
% zchan is optional. If dde is not yet initiated, initiate here.

% check if need to init
if nargin < 2,
   zchan = ddeinit('zemax','topic');
   retstr = ddereq(zchan,'GetRefresh',[1 1]);
end

% get number and names of surfaces
surfs = fieldnames(Lens);
nlensurf = length(surfs);

% get # of surfaces of lens currently in buffer
sysstr = ddereq(zchan,'GetSystem',[1 1]);
ncursurf = sscanf(sysstr,'%d',1); % doesn't count the object surface #0

% insert or delete surfaces as necessary
if ncursurf+1 > nlensurf,
   for n = nlensurf+1:ncursurf,
      retstr = ddereq(zchan,['DeleteSurface,' num2str(n)],[1 1]);
   end
elseif ncursurf+1 < nlensurf,
   for n = ncursurf+1:nlensurf,
      retstr = ddereq(zchan,['InsertSurface,' num2str(n)],[1 1]);
   end
end

% enter each surface
for n = 1:nlensurf,
   eval(['surftmp = Lens.' surfs{n} ';']);
   z_setsurf(zchan,surftmp);
end

% push the lens to the ZEMAX Lens Data Editor
retstr = ddereq(zchan,'PushLens',[1 1],10000);
status = sscanf(retstr,'%d',1);
if status ~= 0, error(['error in PushLens: ' num2str(status)]); end

return


