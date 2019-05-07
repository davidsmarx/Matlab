function [Im, sP] = ReadCodeVExcelImage(fn)
% [Im, sP] = ReadCodeVExcelImage(fn)
%
% fn must be full path
%
% 2019-05-07
% initial hack attempt for reading Excel files provided by Blake
% intend to evolve, especially if I learn CodeV

U = CConstants;

[Excel, Workbook, Sheets, sheet] = Exlopen(fn);

% header
irow = 1;
atmp = '';
while ~strcmp(atmp,'Array Size:') && irow < 22,
    irow = irow + 1;
    atmp = Exlgetval(sheet, ['a' num2str(irow)]);
    if ischar(atmp),
        atmp = strtrim(atmp);
        switch atmp
            case 'Grid spacing:'
                sP.dx = Exlgetval(sheet, ['b' num2str(irow)])*U.MM;
                sP.dy = Exlgetval(sheet, ['d' num2str(irow)])*U.MM;
        end
    end
end

sP.nx = Exlgetval(sheet, ['b' num2str(irow)]);
sP.ny = Exlgetval(sheet, ['c' num2str(irow)]);

% read the data
Im = Exlgetval(sheet, {['b' num2str(irow+1)], [num2column(1+sP.nx) num2str(irow+sP.ny)]});

Exlclose(Excel);

