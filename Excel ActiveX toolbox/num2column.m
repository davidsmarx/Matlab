function strCol = num2column(num)
% column = num2column(num)

if num <= 0, error('column starts with 1'); end
if num == 1, strCol = 'A'; return, end

% number of digits
ndig = ceil(log(num)/log(26));

numrem = num-1;
for ii = 1:ndig,
    bb = (26.^(ndig-ii));
    column(ii) = floor((numrem)/bb);
    numrem = rem((numrem),bb);
    sCol(ii) = 'A' + column(ii);
end

disp(column);

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
