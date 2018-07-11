%   Copyright 2016, 2017 California Institute of Technology
%   Users must agree to abide by the restrictions listed in the
%   file "LegalStuff.txt" in the PROPER library directory.
%
%   PROPER developed at Jet Propulsion Laboratory/California Inst. Technology
%   Original IDL version by John Krist
%   Matlab translation by Gary Gutt


function apm = prop_ellipse(bm, rx, ry, varargin)
%        apm = prop_ellipse(bm, rx, ry, varargin)
% Return an image containing an antialiased, filled ellipse
%
% Outputs:
% apm  = aperture mask containing antialiased filled ellipse
%
% Required Inputs:
% bm   = beam structure
% rx   = radius along x (meters unless norm = 1, then fraction of beam radius)
% ry   = radius along y (meters unless norm = 1, then fraction of beam radius)
%
% Optional inputs:
% 'xc'                = center of ellipse relative to wf center X
%                       (m unless 'norm')
% 'yc'                = center of ellipse relative to wf center Y
%                       (m unless 'norm')
% 'dark'              : draw a dark circle (0 inside, 1 outside)
%                       (default is opposite way)
% 'norm'              = 1 radii and center coordinates are normalized
%                       to beam radius

% 2005 Feb     jek  created idl routine
% 2008 Sep     jek  Fixed bug that was causing a dark line when the
%                   ellipse edge was within 1/75 of the inner pixel boundary.
% 2014 May 13  gmg  Matlab translation
% 2015 Nov 17  gmg  Fixed bug in xmin, xmax to fill in ellipse correctly.
% 2016 Jun 14  gmg  Fixed bug for case where cx is outside of grid
% 2017 Mar 01  gmg  Revised for keyword/value for optional inputs
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Set default values of input parameters
  cx   = 0.0;                           % center of rectangle X
  cy   = 0.0;                           % center of rectangle Y
  dark = 0;
  norm = 0;

% Set values of internal parameters
  del  = 0.0000001;
  dx   = prop_get_sampling(bm);         % spacing between points in x (m)
  dy   = prop_get_sampling(bm);         % spacing between points in y (m)
  [ny, nx] = size(bm.wf);               % number of pixels in wavefront array
  prx  = prop_get_beamradius(bm) / dx;  % beam radius x in pixels
  pry  = prop_get_beamradius(bm) / dy;  % beam radius y in pixels

  icav = 0;                     % index in cell array varargin
  while icav < size(varargin, 2)
    icav = icav + 1;
    switch lower(varargin{icav})
      case {'cx', 'xc'}
        icav = icav + 1;
        cx   = varargin{icav};
      case {'cy', 'yc'}
        icav = icav + 1;
        cy   = varargin{icav};
      case {'dark'}
        dark = 1;
      case {'norm'}
        icav = icav + 1;
        norm = varargin{icav};
      otherwise
        error('prop_ellipse: Unknown keyword: %s\n', varargin{icav});
    end
  end

  cpx  = floor(nx / 2 + 1) + cx / dx;   % center X (pixels)
  cpy  = floor(ny / 2 + 1) + cy / dy;   % center Y (pixels)
  radx = rx / dx;                       % radius in X (pixels)
  rady = ry / dy;                       % radius in Y (pixels)

  if norm == 1
    cpx  = floor(nx / 2 + 1) + cx * prx;
    cpy  = floor(ny / 2 + 1) + cy * pry;
    radx = rx * prx;                    % radius in X (pixels)
    rady = ry * pry;                    % radius in Y (pixels)
  end

  npe  = 1000;                          % number of points in ellipse
  npe1 = npe - 1;
  ang  = 2.0 * pi * [0 : npe1] / npe1;  % angle (rad)
  elx  = radx * cos(ang);               % ellipse coordinates x (pixels)
  ely  = rady * sin(ang);               % ellipse coordinates y (pixels)
  difx = elx(2 : npe) - elx(1 : npe1);  % coordinate differences x (pixels)
  dify = ely(2 : npe) - ely(1 : npe1);  % coordinate differences y (pixels)
  dang = ang(2) / max(sqrt(difx.^2 + dify.^2)) / 100.0;

  npe  = floor(2 * pi / dang);          % number of points in ellipse
  npe1 = npe - 1;
  ang  = 2.0 * pi * [0 : npe1] / npe1;  % angle (rad)
  elx  = radx * cos(ang) + cpx;         % ellipse coordinates x (pixels)
  ely  = rady * sin(ang) + cpy;         % ellipse coordinates y (pixels)
% Allow for ellipse with portions outside of grid
  elx  = max(elx,  1);
  elx  = min(elx, nx);
  ely  = max(ely,  1);
  ely  = min(ely, ny);

  ix   = find(elx - fix(elx) ==  0.5);
  elx(ix)  = elx(ix) - del;
  ix   = find(elx - fix(elx) == -0.5);
  elx(ix)  = elx(ix) + del;

  iy   = find(ely - fix(ely) ==  0.5);
  ely(iy)  = ely(iy) - del;
  iy   = find(ely - fix(ely) == -0.5);
  ely(iy)  = ely(iy) + del;

  elx  = round(elx);
  ely  = round(ely);

  toto = sum(elx >= 1 & elx <= nx & ely >= 1 & ely <= ny);
  if toto == 0
    apm  = ones(ny, nx);
    return
  end

  apm  = zeros(ny, nx);
  apm(sub2ind([ny, nx], ely, elx)) = 1.0;

  nsub = 11;
  nsbn = -floor(nsub / 2);
  nsbp =   ceil(nsub / 2) - 1;
  [x0, y0] = meshgrid([nsbn : nsbp], [nsbn : nsbp]);
  x0   = x0 / nsub;
  y0   = y0 / nsub;

  nsb2 = nsub^2;
  rdx2 = radx^2;
  rdy2 = rady^2;

  for iy = 1 : ny
    for ix = 1 : nx
      if apm(iy, ix) == 1.0
        subx = x0 + ix;
        suby = y0 + iy;
        pix  = (double((subx - cpx).^2 / radx^2) ...
              + double((suby - cpy).^2 / rady^2)) <= 1.0;
        apm(iy, ix) = max((sum(sum(pix)) / nsb2), 1e-8);
      end
    end
  end

% Fill in the aperture area row-by-row
  cpx  = max(1, min(round(cpx), nx));
  for iy = max(min(ely) + 1, 1) : min(max(ely) - 1, ny)
    xmin = max( 1, max(find(apm(iy,  1 : cpx))) + 1);
    xmax = min(nx, min(find(apm(iy, cpx : nx))) + cpx - 2);
    apm(iy, xmin : xmax) = 1.0;
  end

  if dark == 1
    apm  = 1.0 - apm;
  end
end                     % function prop_ellipse
