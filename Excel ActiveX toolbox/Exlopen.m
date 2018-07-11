function [Excel, Workbook, Sheets, Activesheet] = Exlopen(filename)
% [Excel, Workbook, Sheets, Activesheet] = Exlopen(filename)
%
% Open Excel, add workbook, change active worksheet,
% filename must include full path
%
% Get a handle to the active sheet
% Activesheet = Excel.Activesheet;

	
% First open an Excel Server
Excel = actxserver('Excel.Application');

try,

set(Excel, 'Visible', 1);
Workbooks = Excel.Workbooks;

if nargin == 0,
	% Insert a new workbook
	Workbook = invoke(Workbooks, 'Add');
else,
	% open workbook
	if exist(filename) == 2 | exist([filename '.xls']) == 2,
		Workbook = invoke(Workbooks,'Open',filename);
	else,
		warning([filename ' does not exist, creating workbook']);
		Workbook = invoke(Workbooks, 'Add');
		invoke(Workbook,'SaveAs',filename);
	end
end

invoke(Workbook,'Activate');

% get the pointer to Sheets
Sheets = Excel.ActiveWorkBook.Sheets;

% get the handle of the activesheet
Activesheet = Excel.Activesheet;

catch,
	% if something goes wrong, release the server
	warning('error opening workbook, releasing Excel');
	release(Excel);
end
	
	
	
return

% get a handle to the second worksheet and activate it
sheet2 = get(Sheets, 'Item', 2);
invoke(sheet2, 'Activate');

% Get a handle to the active sheet
Activesheet = Excel.Activesheet;

% Put a MATLAB array into Excel
A = [1 2; 3 4];  
ActivesheetRange = get(Activesheet,'Range','A1','B2');
set(ActivesheetRange, 'Value', A);

% Get back a range.  It will be a cell array, since the cell range can
% contain different types of data.
Range = get(Activesheet, 'Range', 'A1', 'B2');
B = Range.value;

% Convert to a double matrix.  The cell array must contain only scalars.
B = reshape([B{:}], size(B));

% Now save the workbook
invoke(Workbook, 'SaveAs', 'myfile.xls');

% To avoid saving the workbook and being prompted to do so,
% uncomment the following code.
% Workbook.Saved = 1;
% invoke(Workbook, 'Close');

% Quit Excel
invoke(Excel, 'Quit');

