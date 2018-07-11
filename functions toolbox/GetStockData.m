function [open, high, low, close, date] = GetStockData(symbol, startdate, frequency, periods)
% [open, high, low, close, date] = GetStockData(symbol, startdate, frequency, periods)
%
% GETSTOCKDATA Fetch historical prices for a given stock symbol
% Inputs:
%   SYMBOL: String representing the stock symbol
%   STARTDATE: String in MM/DD/YYYY format, or serial date number
%   FREQUENCY: Daily ('d'), Monthly ('m'), or Yearly ('y')
%   PERIOD: Number of data points to return (number of days, months or
%           years)
% Outputs:
% The outputs are all parallel column arrays, with PERIODS number of
% entries. Each day 
%   OPEN: Stock opening price
%   HIGH: High price during that trading day
%   LOW: Low price during that trading day
%   CLOSE: Closing market price
%   DATE: The date for which the data is valid

% If the input date is a string, convert it to a serial date number
if (isstr(startdate))
    startdate = datenum(startdate);
end

% Query the server for the stock data. Return the stock data as a 
% cell array of strings.
stockdata = QueryServerForData(symbol, startdate, frequency, periods);

% Make the cell array of strings into one long string. The result is 
% a comma-separated list of values, grouped in sets of six: DATE,
% OPEN, HIGH, LOW, CLOSE, VOLUME. This pattern repeats through the
% string.
stockdata = cat(2, stockdata{:});

% Parse the string data into MATLAB numeric arrays. 
try,
   [dates, open, high, low, close] = ...
      strread(stockdata,'%s%f%f%f%f%*n', 'delimiter', ',', 'emptyvalue', NaN);
catch,
   disp('stockdata return:');
   disp(stockdata);
   keyboard;
end

% Convert the string dates into date numbers. The plotting functions
% need these date numbers.
date = datenum(dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stockdata = QueryServerForData(symbol, startdate, frequency, periods)
% QUERYSERVERFORDATA Return stock data from an Internet server as a cell array
% This function uses the historical stock price server at
% finance.yahoo.com. 
% Inputs:
%   SYMBOL: String representing the stock symbol
%   STARTDATE: Serial date number
%   FREQUENCY: Daily ('d'), Monthly ('m'), or Yearly ('y')
%   PERIOD: Number of data points to return (number of days, months or
%           years)

% Build the query string, a specially formatted URL 
urlString = BuildUrlString(symbol, startdate, frequency, periods);

% Send the query to the server by opening the URL. Use the Java URL
% class to establish the connection.
disp('Contacting server...');
url = java.net.URL(urlString);

% Once the connection is established, create a stream object to read the
% data that the server has returned. Use the stream object to create a
% buffered i/o stream object, which allows the use of readLine -- this
% enables us to read an entire line at a time, rather than just a single
% character.
stream = openStream(url);
ireader = java.io.InputStreamReader(stream);
breader = java.io.BufferedReader(ireader);

% Skip the first line, which consists only of column header labels
line = readLine(breader);

% Read all the available data. We know we've come to the end of the data
% when the readLine call returns a zero length string. Store the data,
% line by line, into a cell array. These strings will eventually be 
% concatenated into one long string, and parsed by STRREAD, so make 
% sure that each string ends with a comma (this uniformity makes for
% easier parsing).

stockdata = {};
disp(['Reading data for symbol ' upper(symbol) '...']);
while 1
    line = readLine(breader);
    if (prod(size(line)) == 0) 
        break;
    end
    line = char(line);
    if line(end)~=',';line(end+1) = ',';end
    % Store each line in a cell array.        
    stockdata{end+1} = line;
end

% Close the streams, in the opposite order in which they were opened.
close(breader);
close(ireader);
close(stream);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = BuildUrlString(symbol, startdate, frequency, periods)
% BUILDURLSTRING Construct a specially formatted URL query string
% Inputs:
%   SYMBOL: String representing the stock symbol
%   STARTDATE: Serial date number
%   FREQUENCY: Daily ('d'), Monthly ('m'), or Yearly ('y')
%   PERIOD: Number of data points to return (number of days, months or
%           years)

% Extract the day, month and year from the starting date
day = datevec(startdate);
startDay = day(3);
startMonth = day(2)-1;
startYear = day(1);
% startDay = day(startdate);
% startMonth = month(startdate);
% startYear = year(startdate);

% Initially, the end date is the same as the starting date
endDay = startDay;
endMonth = startMonth;
endYear = startYear;

% Given the frequency (daily, monthly, yearly), compute the end date
switch frequency
    case 'y' 
        endYear = endYear + periods;
    case 'm'
        endMonth = endMonth + (rem(periods, 12));
        endYear = endYear + (floor(periods./12));
    case 'd'
        newdate = startdate + periods;
        newday = datevec(newdate);
        endDay = newday(3);
        endMonth = newday(2)-1;
        endYear = newday(1);
%         endDay = day(newdate);
%         endMonth = month(newdate);
%         endYear = year(newdate);

end        

% The server name is constant, but the query data tagged onto the end of it
% varies according to user input
server = 'http://table.finance.yahoo.com/table.csv';
query = ['?a=' num2str(startMonth) '&b=' num2str(startDay) '&c=' num2str(startYear) ...
         '&d=' num2str(endMonth)   '&e=' num2str(endDay)   '&f=' num2str(endYear) ...
         '&s=' num2str(symbol) '&y=0&g=' num2str(frequency) ];

% Concatenate the server name and query data together into the complete
% query string
s = [ server query ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
