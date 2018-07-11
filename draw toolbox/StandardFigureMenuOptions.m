function StandardFigureMenuOptions(hfig)

% file menu
% first delete the built in file menu
hmnufile = findall(hfig, 'tag', 'figMenuFile');
delete(hmnufile);
hmnufile = uimenu(hfig, 'Label', 'File', 'Tag', 'mnuFile');

% allow copy to clipboard (from Matlab Technical Solutions)
h = uimenu(hfig, 'Label', 'Edit');
uimenu(h, 'Label', 'Copy Figure', 'Callback', {@localprint});

% % allow zooming menu options (from Matlab Technical Solutions)
% f = uimenu(hfig, 'Label', 'Tools');
% uimenu(f,'Label','Zoom In','Callback','zoom(2)', 'Accelerator', 'z');
% uimenu(f,'Label','Zoom Out','Callback','zoom(1/2)', 'Accelerator', 'o');

function localprint(gcbo, eventdata)
%print(gcf, '-dmeta');
set(gcf, 'InvertHardcopy', 'off');
print(gcf, '-dbitmap');

