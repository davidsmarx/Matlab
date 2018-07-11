function y = decibel(x)
% y = decibel(x)
% y = 10.*log10(abs(x))

if isa(x,'sym'),
   y = 10*log10(abs(x));
else,
ax = abs(x);
nn = ax>1e-10;
y = zeros(size(ax));
y(nn) = 10.*log10(ax(nn));
y(~nn) = -100;
end

return