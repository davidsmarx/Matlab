

b = (peaks+400)*10^(-30);  % A sample surface data with very small range, riding on top of a large value

figure; surf(b) % No surface shown!

set(gcf, 'Renderer', 'zbuffer') % Surface shown!


% %% Alternative : Normalizing data
% lowest = min(min(b));  % Calculate the lowest point
% zeroed_b = (b - lowest); % Shift the surface so that lowest point is zero
% normalized_b = zeroed_b / (max(max(zeroed_b))); % Normalize so that the highest point is 1
% 
% figure; surf(normalized_b)  % Some surface shown

