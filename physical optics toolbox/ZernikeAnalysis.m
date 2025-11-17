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
% CheckOption('x', [], varargin{:});
% CheckOption('y', [], varargin{:});
% CheckOption('com', [0 0], varargin{:});
% CheckOption('Rnorm', [], varargin{:});
% CheckOption('modes', 1:36, varargin{:});
% CheckOption('polyorder', 'Noll', varargin{:});
% CheckOption('do_phaseunwrap', true
%
% return:
% ZZ = Zernike coeffs, (rad normalized rms)
% phaimg = phaimg(bMask) is input to zernikeval and then ptt removed
% phares = residual phase = bMask .* (phaimg - zernikeval(ZZ))
% sOptions = struct(
%    'bMask', bMask,
%    'phaimg_input', phaimg_input, = the unadulterated phase map for Zernike fitting
%    'Rnorm', Rnorm,
%    'Nz', Nz,
%    'xim', xim, 'yim', yim,
%    'phafit', phafit, = Zernike fit including ptt, no bMask
%    'phafit_ptt', phafit_ptt, = Zernike fit without ptt, no bMask


% options
bMask = CheckOption('bMask', [], varargin{:});
bField_is_phase = CheckOption('isphase', false, varargin{:}); % if field = real-valued phase array
Rnorm = CheckOption('Rnorm', [], varargin{:});
Nz = CheckOption('Nz', 1:36, varargin{:}); % either 'Nz' or 'modes'
Nz = CheckOption('modes', Nz, varargin{:});
polyorder = CheckOption('polyorder', 'Noll', varargin{:});
do_phaseunwrap = CheckOption('do_phaseunwrap', true, varargin{:});
xim = CheckOption('x', [], varargin{:});
yim = CheckOption('y', [], varargin{:});
com = CheckOption('com', [0 0], varargin{:});

% if input is phase, must include bMask
if bField_is_phase && isempty(bMask)
    error('if input is phase, must include bMask');
end

%
if isempty(bMask),
    [~, bMask] = AutoMetric(abs(field));
end

% make a grid
if isempty(xim) || isempty(yim)
    [xim, yim] = CreateGrid(bMask);
end
xim = xim - com(1);
yim = yim - com(2);
[Xim, Yim] = meshgrid(xim, yim);
Rim = hypot(Xim, Yim);

% phase
if bField_is_phase
    phaimg = field;
else
    phaimg = angle(field);
end

phaimg_input = phaimg; % save the input phase as it was input

if do_phaseunwrap
    phaimg = unwrap_HCIT(phaimg, bMask, 'selem', []); % don't alter the bMask
end

% remove piston (no mod2pi() because it is unwrapped phase) and apply bMask
phaimg(~bMask) = 0;
phaimg(bMask) = phaimg(bMask) - mean(phaimg(bMask));

% norm radius for zernikefit
if isempty(Rnorm),
    Rnorm = max(Rim(bMask(:)));
end

ZZ = zernikefit(Xim(bMask), Yim(bMask), phaimg(bMask), Nz, Rnorm, polyorder);
ZZ_ptt = ZZ; ZZ_ptt(1:3) = 0;

% the fit:
phafit = zernikeval(ZZ, Xim, Yim, Rnorm, polyorder);
phafit_ptt = zernikeval(ZZ_ptt, Xim, Yim, Rnorm, polyorder);

% residuals of the Zernike fit
phares = zeros(size(phaimg));
phares(bMask) = phaimg(bMask) - zernikeval(ZZ, Xim(bMask), Yim(bMask), Rnorm, polyorder);

% remove ptt from phaimg
phaimg(bMask) = phaimg(bMask) - zernikeval(ZZ(1:3), Xim(bMask), Yim(bMask), Rnorm, polyorder);
phaimg(~bMask) = 0; % redundant

% return values:
% ZZ, phaimg, phares, sOptions
sOptions = struct('bMask', bMask, 'phaimg_input', phaimg_input, ...
    'Rnorm', Rnorm, 'Nz', Nz,...
    'xim', xim, 'yim', yim, ...
    'phafit', phafit, 'phafit_ptt', phafit_ptt);

end % main

function phunwrap = unwrap_(pha, bMask)

phunwrap = pha;

end

