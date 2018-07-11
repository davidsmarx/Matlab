function cTickLabel = XTickLabelPi(hax)
% cTickLabel = XTickLabelPi(hax)
% append string '\pi' to x tick labels

cTickLabel = get(hax,'XTickLabel');
% cell array of tick label strings

for il = 1:length(cTickLabel),
    cTickLabel{il} = [cTickLabel{il} '\pi'];
end

set(gca,'XTickLabel',cTickLabel)

