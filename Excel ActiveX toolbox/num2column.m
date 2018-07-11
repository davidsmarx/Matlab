function column = num2column(num)
% column = num2column(num)

if num <= 0, error('column starts with 1'); end

first = floor((num-1)/26);
sec   = rem((num-1),26) + 1;

if first < 1,
    column = 'A' + sec - 1;
else,
    column = [ 'A' + first - 1, 'A' + sec - 1 ];
end

column = char(column);
