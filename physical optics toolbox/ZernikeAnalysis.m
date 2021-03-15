function [ZZ, phaimg, phares, sOptions] = ZernikeAnalysis(field, varargin)
% [ZZ, phaimg, phares, sOptions] = ZernikeAnalysis(field, varargin)
%
% wrapper for zernikefit
%
% field = complex array. We apply zernike fit to phase
%
% options:
% CheckOption('bMask', [], varargin{:});
% CheckOption('Rnorm', [], varargin{:});
% CheckOption('modes', 1:36, varargin{:});
% CheckOption('polyorder', 'Noll', varargin{:});
%
% return:
% ZZ = Zernike coeffs, (rad normalized rms), ZZ(1:3) == 0
% phaimg = phase with ZZ(1:3) = 0, (rad), phaimg(~bMask) = 0
% phares = residual phase = phaimg - zernikeval(ZZ)
% sOptions = struct('bMask', bMask, 'Rnorm', Rnorm, 'Nz', Nz);


% options
bMask = CheckOption('bMask', [], varargin{:});
Rnorm = CheckOption('Rnorm', [], varargin{:});
Nz = CheckOption('Nz', 1:36, varargin{:}); % either 'Nz' or 'modes'
Nz = CheckOption('modes', Nz, varargin{:});
polyorder = CheckOption('polyorder', 'Noll', varargin{:});

%
if isempty(bMask),
    [~, bMask] = AutoMetric(abs(field));
end

% make a grid
[xim, yim, Xim, Yim, Rim] = CreateGrid(bMask);

% phase and adjust piston
phaimg = angle(field);
phaimg(bMask) = mod2pi(phaimg(bMask) - mean(phaimg(bMask)));

% norm radius for zernikefit
if isempty(Rnorm),
    Rnorm = max(Rim(bMask));
end

ZZ = zernikefit(Xim(bMask), Yim(bMask), phaimg(bMask), Nz, Rnorm, polyorder);

% remove ptt from phaimg
phaimg(bMask) = phaimg(bMask) - zernikeval(ZZ(1:3), Xim(bMask), Yim(bMask), Rnorm, polyorder);
phaimg(~bMask) = 0;

% residual of the rest of the moes
phares = zeros(size(phaimg));
phares(bMask) = phaimg(bMask) - zernikeval([zeros(3,1); ZZ(4:end)], Xim(bMask), Yim(bMask), Rnorm, polyorder);

% return values:
% ZZ, phaimg, phares, sOptions
sOptions = struct('bMask', bMask, 'Rnorm', Rnorm, 'Nz', Nz);

