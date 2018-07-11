function make_figure(H,fn, FSize)
%return
if ~exist('FSize')
    FSize = [0 0 6 4];
end
    set(H,'PaperUnits','inches','PaperPosition',FSize);
    %exportfig(H, [fn '.ps'], 'Format', 'eps', 'Color', 'cmyk','FontSize', 1.2);
    %exportfig(H, [fn '.ps'], 'Format', 'eps2', 'Color', 'cmyk','FontSize', 1.2);
    %exportfig(H, [fn '.ps'], 'Format', 'eps2', 'Resolution',30,'Color', 'cmyk','FontSize', 1.2);
    %exportfig(gcf,fn,'Resolution',100);
    print('-dpng',[fn '.png']);  
    crop_fig([fn '.png']);
    %print('-dpng',[fn '.png'], '-r200');  
    %system(sprintf('/home/bseo/utils/ps2eps -f -F %s.ps',fn));
    %system(sprintf('/usr/local/bin/epstopdf %s.eps',fn));
    %system(sprintf('rm -f %s.ps',fn));
    %system(sprintf('rm -f %s.eps',fn));
return

