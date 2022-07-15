function save_as_png(handle, filename, r1)
% save_as_png(handle, filename, r1)
%   r1 = 300 (default) sets '-r300' sets the output resolution to 300 dots per inch
%        0 specifies screen resolution, '-r0'. 

% default resolution
if ~exist('r1','var'),
    r1 = 300; 
end

    lbwh = get(handle,'Position');
    r0=get(0,'ScreenPixelsPerInch');
    pos = [0 0 lbwh(3) lbwh(4)]/r0;
    
    set(handle, 'PaperUnits', 'inches','PaperSize', pos(3:4),'PaperPosition', pos);
    print(handle, '-dpng', ['-r' num2str(r1)], filename);
end