function ExlCopyFigureToSheet(hfig, sheet, cellid)
% ExlCopyFigureToSheet(hfig, sheet, cellid)

figure(hfig); drawnow;
set(hfig,'invertHardcopy','off')
print -dbitmap

rh = sheet.Range(cellid);
sheet.Paste(rh);

