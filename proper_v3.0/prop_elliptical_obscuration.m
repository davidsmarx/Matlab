%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function bm = prop_elliptical_obscuration(bm, rx, ry, varargin)
%        bm = prop_elliptical_obscuration(bm, rx, ry, varargin)
% Multiplies the wavefront in bm by an elliptical obscuration
%
% Outputs:
% bm   = beam structure (output)
%
% Required inputs:
% bm   = beam structure (input)
% rx   = radius of obscuration X (meters unless 'norm' is set)
% ry   = radius of obscuration Y (meters unless 'norm' is set)
%
% Optional inputs:
% 'xc'                = center coordinate X (meters unless 'norm')
% 'yc'                = center coordinate Y (meters unless 'norm')
% 'norm'              : for rx, ry, xc, yc normalized to beam radius

% 2005 Feb     jek  created idl routine
% 2016 Apr 14  gmg  Matlab translation
% 2017 Mar 01  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  cx   = 0.0;                   % obscuration center to beam center X
  cy   = 0.0;                   % obscuration center to beam center Y
  norm = 0;

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};  % obscuration center to beam center X
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};  % obscuration center to beam center Y
      case {'norm'}
        norm = 1;
      otherwise
        error('prop_elliptical_obscuration: Unknown keyword: %s\n', ...
          varargin{icav});
    end
  end

  bm.wf = bm.wf .* prop_shift_center(prop_ellipse(bm, rx, ry, ...
          'cx', cx, 'cy', cy, 'norm', norm, 'dark'), 'inv');
end                     % function prop_elliptical_obscuration
