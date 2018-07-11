function mask = CreateMaskFromAmp(amp)
% mask = CreateMaskFromAmp(amp)

[cnt, x_amp] = hist(amp(:),21);

is_min = [false cnt(2:end-1) <= cnt(1:end-2) & cnt(2:end-1) <= cnt(3:end) false];
if ~any(is_min),
    error('no threshold found');
end

minima = x_amp(is_min);
thresh = minima(1);

mask = amp >= thresh;