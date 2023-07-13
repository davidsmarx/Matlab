function [xfield, yfield, weight, vdx, vdy, vcx, vcy, van] =...
   z_getfield(zchan,nfield)
% [xfield, yfield, weight, vdx, vdy, vcx, vcy, van] =...
%   z_getfield(zchan,nfield)
%
% nfield = a valid field number
% [xfield, yfield, ...] correspond to the "Field Data" form in ZEMAX
%
% [type, nn, field_max_xy, norm_method] = z_getfield(zchan)
%
% type: 0 => angle in rad (Zemax returns deg),
%       1 => object height in m
%       2 => paraxial image
%       3 => real image height
%       4 => theodolite angles
%
% nn = number of field points currently defined
% field_max_xy = [x y] = normalization radius used to calculate [hx hy]
% norm_method = 'radial' or 'rectangular'

unitsdefinitions;

% first check number and type of field points
retstr = ddereq(zchan,'GetField, 0',[1 1]);
if ~ischar(retstr),
    % something is wrong
    [xfield, yfield, weight, vdx, vdy, vcx, vcy, van] = deal([]);
    return
end

retval = sscanf(retstr, '%f,',[1 inf]);
type = retval(1); % type
nn   = retval(2); % nn

switch nargin,

    case 1,
        % make recursive call to get the field points and determine the
        % maximum field radius

        retstr = ddereq(zchan,'GetField, 0',[1 1]);
        cc = textscan(retstr,'%d,%d,%f,%f,%d');
        [type, number, max_x, max_y, norm_method] = deal(cc{:});
        switch norm_method
            case 0
                norm_method = 'radial';
            case 1
                norm_method = 'rectangular';
        end
        % the outputs:
        [xfield, yfield, weight, vdx] = ...
            deal(type, number, ApplyUnits([max_x max_y], type), norm_method);
        
    case 2,

        if nfield <= 0 | nfield > nn,
            error('invalid nfield');
        end

        % send the DDE command
        cmdstr = sprintf('GetField,%d',nfield);
        retval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

        % parse results
        xfield = ApplyUnits(retval(1),type);
        yfield = ApplyUnits(retval(2),type);
        weight = retval(3);
        vdx = retval(4);
        vdy = retval(5);
        vcx = retval(6);
        vcy = retval(7);
        van = retval(8);

    otherwise,
        disp(['usage: [xfield, yfield, weight, vdx, vdy, vcx, vcy, van] ='...
            'z_getfield(zchan,nfield)']);

end
        
function fieldmks = ApplyUnits(fieldin,type)

unitsdefinitions;

switch type,
    case 0,
        fieldmks = fieldin*P;
    case 1,
        fieldmks = fieldin*MM;
    case 2,
        fieldmks = fieldin*MM;
    otherwise,
        error('uknown type');
end
