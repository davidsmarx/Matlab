function strCol = num2column(num)
% column = num2column(num)

if num <= 0, error('column starts with 1'); end
%if num == 1, strCol = 'A'; return, end

% number of digits
ndig = ceil(log(num)/log(26));

numrem = num;

% ones-place
column = rem(numrem - 1, 26) + 1;
sCol = 'A' + column - 1;
numrem = floor((numrem - column) / 26);

% higher digits
while numrem > 0,
    thiscolumn = rem(numrem - 1, 26) + 1;
    column = [thiscolumn column];
    sCol = ['A' + thiscolumn - 1, sCol];
    numrem = floor((numrem - thiscolumn) / 26);
end

%disp(column);

% 
% first = floor((num-1)/26);
% sec   = rem((num-1),26) + 1;
% 
% if first < 1,
%     column = 'A' + sec - 1;
% else,
%     column = [ 'A' + first - 1, 'A' + sec - 1 ];
% end
% 
strCol = char(sCol);
