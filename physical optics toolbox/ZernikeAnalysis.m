function [ZZ, phaimg, phares, sOptions] = ZernikeAnalysis(field, varargin)
% [ZZ, phaimg, phares, sOptions] = ZernikeAnalysis(field, varargin)
%
% wrapper for zernikefit
%
% field = complex array. We apply zernike fit to phase
%
% options:
% CheckOption('bMask', [], varargin{:}); % default: [~, bMask] = AutoMetric(abs(field));
% CheckOption('isphase', false, varargin{:}); % if field = real-valued phase array
% CheckOption('Rnorm', [], varargin{:});
% CheckOption('modes', 1:36, varargin{:});
% CheckOption('polyorder', 'Noll', varargin{:});
% CheckOption('do_phaseunwrap', true
%
% return:
% ZZ = Zernike coeffs, (rad normalized rms), ZZ(1:3) == 0
% phaimg = phase with ZZ(1:3) = 0, (rad), phaimg(~bMask) = 0
% phares = residual phase = phaimg - zernikeval(ZZ)
% sOptions = struct('bMask', bMask, 'Rnorm', Rnorm, 'Nz', Nz);


% options
bMask = CheckOption('bMask', [], varargin{:});
bField_is_phase = CheckOption('isphase', false, varargin{:}); % if field = real-valued phase array
Rnorm = CheckOption('Rnorm', [], varargin{:});
Nz = CheckOption('Nz', 1:36, varargin{:}); % either 'Nz' or 'modes'
Nz = CheckOption('modes', Nz, varargin{:});
polyorder = CheckOption('polyorder', 'Noll', varargin{:});
do_phaseunwrap = CheckOption('do_phaseunwrap', true, varargin{:});

% if input is phase, must include bMask
if bField_is_phase && isempty(bMask)
    error('if input is phase, must include bMask');
end

%
if isempty(bMask),
    [~, bMask] = AutoMetric(abs(field));
end

% make a grid
[xim, yim, Xim, Yim, Rim] = CreateGrid(bMask);

% phase
if bField_is_phase
    phaimg = field;
else
    phaimg = angle(field);
end

if do_phaseunwrap
    phaimg = unwrap_HCIT(phaimg, bMask, 'selem', []); % don't alter the bMask
end

% remove piston (no mod2pi() because it is unwrapped phase)
phaimg(bMask) = phaimg(bMask) - mean(phaimg(bMask));

% norm radius for zernikefit
if isempty(Rnorm),
    Rnorm = max(Rim(bMask(:)));
end

ZZ = zernikefit(Xim(bMask), Yim(bMask), phaimg(bMask), Nz, Rnorm, polyorder);

% the fit:
phafit = zernikeval(ZZ, Xim, Yim, Rnorm, polyorder);

% remove ptt from phaimg
phaimg(bMask) = phaimg(bMask) - zernikeval(ZZ(1:3), Xim(bMask), Yim(bMask), Rnorm, polyorder);
phaimg(~bMask) = 0;

% residual of the rest of the moes
phares = zeros(size(phaimg));
phares(bMask) = phaimg(bMask) - zernikeval([zeros(3,1); ZZ(4:end)], Xim(bMask), Yim(bMask), Rnorm, polyorder);

% return values:
% ZZ, phaimg, phares, sOptions
sOptions = struct('bMask', bMask, 'Rnorm', Rnorm, 'Nz', Nz, 'xim', xim, 'yim', yim, 'phafit', phafit);

end % main

function phunwrap = unwrap_(pha, bMask)

phunwrap = pha;

end

