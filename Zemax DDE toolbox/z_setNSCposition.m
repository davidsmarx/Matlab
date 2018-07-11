function opval = z_setNSCposition(zchan, nsurf, nobject, varargin)
% opval = z_setNSCparm(zchan, nsurf, nobject, position, material)
% opval = z_setNSCparm(zchan, nsurf, nobject, struct_position)
%
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC
% position = [x, y, z, x_tilt, y_tilt, z_tilt], or
% struct_position = struct with fields:
%    x, y, z, x_tilt, y_tilt, z_tilt, material

UU = CConstants;

switch length(varargin),
    case 0,
        disp('zero!');
        
    case 1,
        if isstruct(varargin{1}),
            x = varargin{1}.x;
            y = varargin{1}.y;
            z = varargin{1}.z;
            x_tilt = varargin{1}.x_tilt;
            y_tilt = varargin{1}.y_tilt;
            z_tilt = varargin{1}.z_tilt;
            material = varargin{1}.material;

        else
            % array of doubles
            [x, y, z, x_tilt, y_tilt, z_tilt] =...
                deal(varargin{1}(1), varargin{1}(2), varargin{1}(3),...
                varargin{1}(4), varargin{1}(5), varargin{1}(6));
            material = '';
        end

    case 2,
        [x, y, z, x_tilt, y_tilt, z_tilt] =...
            deal(varargin{1}(1), varargin{1}(2), varargin{1}(3),...
            varargin{1}(4), varargin{1}(5), varargin{1}(6));
        material = varargin{2};
        
    otherwise,
        error('usage');
end

cmdstr = ['SetNSCPosition,' num2str(nsurf) ',' num2str(nobject) ','];
SendValue([cmdstr '1,' num2str(x/UU.MM,'%e')]);
SendValue([cmdstr '2,' num2str(y/UU.MM,'%e')]);
SendValue([cmdstr '3,' num2str(z/UU.MM,'%e')]);
SendValue([cmdstr '4,' num2str(x_tilt/UU.P,'%e')]);
SendValue([cmdstr '5,' num2str(y_tilt/UU.P,'%e')]);
SendValue([cmdstr '6,' num2str(z_tilt/UU.P,'%e')]);
SendValue([cmdstr '7,' material]);

    function SendValue(cmdstr)
        opval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);
    end

end % main
