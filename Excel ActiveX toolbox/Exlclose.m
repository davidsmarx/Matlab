function Exlclose(Excel, varargin)
% Exlclose(Excel, sOptions)
%
% sOptions = struct(...
%     'Alerts', (default = true)
%     )

sOptions = ValidateOptions(varargin{:});

try,
    if ~sOptions.Alerts,
         Excel.EnableEvents = false;
         Excel.DisplayAlerts = false;
    end
    
   % Quit Excel
   Excel.Quit;

   % delete the process
   delete(Excel);
   
catch,
   warning('error produced attempting to quit Excel');
   
end

return

function sOptions = ValidateOptions(varargin)

% default options
sOptions = struct(...
    'Alerts', true ...
    );

% if options are specified:
if nargin >= 1,
    sIn = varargin{1};
    if isstruct(sIn),
        fnames = fieldnames(sIn);
        for ii = 1:length(fnames),
            sOptions.(fnames{ii}) = sIn.(fnames{ii});
        end
    end
end
