function retval = z_setfield(zchan, varargin)
% status = z_setfield(zchan, field_struct)
% status = z_setfield(zchan, n, xf, yf, wgt, vdx, vdy, vcx, vcy, van)
% status = z_setfield(zchan, 0, type, number, normalization)
%
% field_struct has fields n, xf, yf, wgt, vdx, vdy, vcx, vcy, van
%                      or n=0, type, number, normalization
%
% xf, yf in [m]
%
% If the value for n is zero, then the field type, total number of fields,
% and field normalization type is set to the new values. The field type is
% 0 for angle, 1 for object height, 2 for paraxial image height, and 3 for
% real image height. The field normalization type is 0 for radial and 1 for
% rectangular.

global MM;

switch length(varargin)
    case 0,
        error('usage');
    case 1,
        if ~isstruct(varargin{1}),
            error('usage');
        end
        sField = varargin{1};
        if sField.n == 0,
            [n, type, number, normalization] = deal(sField.n, sField.type, sField.number, sField.normalization);
        else
            [n, xf, yf, wgt, vdx, vdy, vcx, vcy, van] = deal(...
                sField.n, sField.xf, sField.yf, sField.wgt,...
                sField.vdx, sField.vdy, sField.vcx, sField.vcy, sField.van);
        end
    otherwise,
        if varargin{1} == 0,
            [n, type, number, normalization] = deal(varargin{:});
        else
            [n, xf, yf, wgt, vdx, vdy, vcx, vcy, van] = deal(varargin{:});
        end
end

if n == 0,
    cmdstr = sprintf('SetField, 0, %d, %d, %d',type,number,normalization);
else
    cmdstr = sprintf('SetField, %d, %f, %f, %f, %f, %f, %f, %f, %f',...
        n, xf/MM, yf/MM, wgt, vdx, vdy, vcx, vcy, van);
end

retval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);

    
% 
% if nargin == 0, disp('usage: status = z_setwave(zchan,lambda)'); return, end
% 
% if ~exist('nwave','var') | isempty(nwave),
%     nwave = 1;
% end
% if ~exist('weight','var') | isempty(weight),
%     weight = 1.0;
% end
% 
% % ZEMAX uses um for wavelength
% cmdstr = sprintf('SetWave,%d,%f,%f',nwave,lam/UM,weight);
% retval = sscanf(ddereq(zchan,cmdstr,[1 1]),'%f,',[1 inf]);
% 
% % call update
% if z_getupdate(zchan), error('z_setwave: GetUpdate Failed!'); end
% 
% return
