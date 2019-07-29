cls

% get result of python model0
pn = '/home/dmarx/HCIT/DST/hcim_model0_run012/results/';
%fn = 'Epups_dmgain_test_20190728T1929.mat';
%fn = 'Epups_dmgain_test_20190728T2007.mat'; % fixed rad-surf, infperact = 8
%fn = 'Epups_dmgain_test_20190728T2115.mat'; % fixed rad-surf, infperact = 7
%fn = 'Epups_dmgain_test_20190728T2127.mat'; % fixed rad-surf, infperact = 6
%fn = 'Epups_dmgain_test_20190728T2150.mat'; % fixed rad-surf, infperact = 6.71
fn = 'Epups_dmgain_test_20190729T0022.mat'; % fixed rad-surf, inf fun scaled by 7.3/10


if ~exist(['./' fn],'file'),
    disp(['copying ' [pn fn] ' to ./']);
    copyfile(PathTranslator([pn fn]), '.');
    disp('done copying');
end

S = load(fn);

dmpk_rad = reshape(S.dmpk_rad, [50 50]);

iwv = 3;
Epupdm1pk = squeeze(S.Epup_dm1pk(iwv,:,:));
Epupdm2pk = squeeze(S.Epup_dm2pk(iwv,:,:));
Epupref   = squeeze(S.Epup_ref(iwv,:,:));

figure_mxn(1,2)
hax(1) = subplot(1,2,1);
imageschcit(angle(Epupdm1pk .* conj(Epupref)))
title('DM1 Poke - Ref')

hax(2) = subplot(1,2,2);
imageschcit(angle(Epupdm2pk .* conj(Epupref)))
title('DM2 Poke - Ref')

set(hax,'clim',[0 1.1])
set(hax,'xlim',150*[-1 1],'ylim',150*[-1 1])
