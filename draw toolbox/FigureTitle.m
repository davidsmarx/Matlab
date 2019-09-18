function han = FigureTitle(stitle, varargin)
% han = FigureTitle(stitle, varargin)
%
% add an annotation text at the top of the figure, 
% useful for adding a single main title to a figure with subplots
%
% varargin can be property, value pairs sent to the annotation handle

han = annotation('textbox', [0.5 0.8 0.2 0.2], 'String', stitle, ...
    'FitBoxToText', 'on', 'LineStyle', 'none', ...
    'FontSize', 24, 'Color', 'r', 'FontWeight', 'bold');
set(han,'HorizontalAlignment','center')
% center horizontally
ppp = get(han,'Position');
set(han,'Position',[0.5 - 0.5*ppp(3) ppp(2:end)])
% so it can be found and deleted later
set(get(han,'parent'),'HandleVisibility','on')

if ~isempty(varargin),
    set(han, varargin{:})
end

