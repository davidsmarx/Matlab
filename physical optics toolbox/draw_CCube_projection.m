function h = draw_CCube_projection(haxes,cc,varargin)
% h = draw_CCube_projection(haxes,cc,options)
% h = draw_CCube_projection(haxes,cc,'lines',options)
%
% cc is a struct of corner cube parameters, see IntMetCom_CornerCube for
% struct definition.
% The fields in the Ccube structure are:
%    size
%    shape ('square' (default), 'triangle', 'sim')
%    edgelen(1:3) = relative length of each roofline. (default = [1 1 1]).
%         The length of each roofline will be size*edgelen(ii).
%    xc  = vertex offset in x (default = 0)
%    yc  = vertex offset in y (default = 0)
%    spin = clocking of corner cube (default = 0)
%    gapwidth(1:3) = size of gap along each roofline (default = 0)
%    dihedral(1:3) = dihedral (rad) error of each roofline
%    rotmatrix(1:3,1:3) = specifies tip/tilt (default = identity matrix)
%    edgedir(1:3,1:3) = matrix of three column vectors, where each vector is a
%       direction vector indicating the direction of a roof line.
%
% if the first option is 'lines', then the corner cube is drawn as lines
% for the perimeter and each roofline. Otherwise, it is drawn using fill3.
% The remainder options are used in the drawing routine.

global MM;

[vp{1}, vp{2}, vp{3}] = CCube2Polygons(cc);

axes(haxes);
C = ['b','r','g'];

if length(varargin) >= 1 & isequal(varargin{1},'lines'),
    for ii = 1:3,
        h{ii} = line(vp{ii}(:,1)/MM,vp{ii}(:,2)/MM,vp{ii}(:,3)/MM,varargin{2:end});
        %h{ii} = patch(vp{ii}(:,1)/MM,vp{ii}(:,2)/MM,vp{ii}(:,3)/MM,C(ii),varargin{:});
    end

else, % use fill3 to draw patches
    for ii = 1:3, % each patch
        for jj = 1:3, % x,y,z
            vpp{jj,ii} = vp{ii}(:,jj)/MM;
        end
        vpp{4,ii} = C(ii);
    end
    h = fill3(vpp{:},varargin{:});
end

view(2); % normal 2-d view
