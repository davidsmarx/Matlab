cls,

filename = 'C:\Tamar\HART\HRT Applications list.xls';

% open an Excel file
Excel = actxserver('Excel.Application');
Excel.Visible = 1;
Workbooks = Excel.Workbooks;
Workbook = Workbooks.Open(filename);
Sheets = Workbook.Sheets;

% the number of sheets in this workbook
nsheets = Sheets.Count;

% get the name of each sheet
for ii = 1:nsheets,
    sheet = Sheets.Item(ii);
    sheetname{ii} = sheet.Name;
end

% user select which sheet to use
disp(strvcat(sheetname{:}));
selectsheet = input('Which Sheet to Activate? ');

% make that sheet active
sheet = Sheets.Item(selectsheet);
sheet.Activate;

% write a value to the active sheet
set(get(sheet,'range','A1'),'value',3.14159);

