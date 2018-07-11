zchan = ddeinit('zemax','topic');
retstr = ddereq(zchan,'GetRefresh',[1 1]);

%pp = 'H:\matlab';
%pp = 'C:\ZEMAX';
%cmdstr = ['GetMetaFile, "' pp '\testmeta.emf", L3d, "' pp '\settings.cfg", 0']
cmdstr = ['GetMetaFile, "H:\ZEMAX lenses\EXPORT.EMF", L3d'];
retstr = ddereq(zchan,cmdstr,[1 1]);
% Lsm   solid model
% Lsh   shaded model
%retstr = ddereq(zchan,'OpenWindow, Lsh',[1 1]);

rc = ddeterm(zchan);
