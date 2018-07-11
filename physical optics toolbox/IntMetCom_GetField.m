function [field, dx, dy, xs, ys, curv] = IntMetCom_GetField(h,bidtmp,varargin)
% [field, dx, dy, xs, ys, curv] = IntMetCom_GetField(h)
% [field, dx, dy, xs, ys, curv] = IntMetCom_GetField(h,BeamID)
% [field, dx, dy, xs, ys, curv] = IntMetCom_GetField(h,BeamID,options)
% [field, dx, dy, xs, ys, curv] = IntMetCom_GetField(h,[],options)
% S = IntMetCom_GetField(...)
%
% IntMetCom_GetField(h,BeamID,options) with no output arguments creates a figure
% with the image of the ampitude of the current field. options are:
%   'dB' or 'abs' (default)
%
% options:
%    'applycurv' (default) or 'noapplycurv'

BeamID = 0;
if exist('bidtmp','var'),
    if ~isempty(bidtmp),
        BeamID = bidtmp;
    end
end
           
[nx, ny, dx, dy, curv, wavelength] = h.GetWavefrontParms(BeamID);

if nx == 0 | ny == 0,
    % wavefront not allocated or BeamID error
    field = 0;
    xs = 0; ys = 0;
    return
end

ptmp = zeros(ny,nx);
[c, d] = h.wavefront_get(BeamID,ptmp,ptmp);
field = c + j*d;

xs = h.x_vector_get(BeamID,zeros(nx,1));
ys = h.y_vector_get(BeamID,zeros(ny,1));

if abs(curv) > 1e-12,
    %warning('non-zero curvature returned');
    if isempty(strmatch('noapplycurv',varargin)),
        [X, Y] = meshgrid(xs,ys);
        R2 = X.^2 + Y.^2; clear X Y;
        atmp = pi*curv/wavelength;
        A = exp(j*atmp*R2);
        clear R2;
        field = A .* field;
        clear A;
        curv = 0;
    end % if applycurv
end

switch nargout,
    case 0,
        MM = 1e-3;

        % parse options
        if strmatch('dB',varargin),
            plfun = @dBa;
            pllim = [-90 0];
        else
            plfun = @abs;
            pllim = [0 max(abs(field(:)))];
        end
        figure, imagesc(xs/MM,ys/MM,plfun(field)), axis image,...
            setimageclim(gcf,pllim), colorbar
        
    case 1,
         S.field = field;
         S.dx = dx;
         S.dy = dy;
         S.xs = xs;
         S.ys = ys;
         S.curv = curv;
         S.wavelength = wavelength;
         
         field = S;
         
    otherwise,
        % normal
end % switch nargout

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

% dx = h.dx;
% dy = h.dy;
% 
% xs = invoke(h,'x_vector_get',zeros(h.nx,1));
% ys = invoke(h,'y_vector_get',zeros(h.ny,1));
% 
% curv = h.curvature;
% 
