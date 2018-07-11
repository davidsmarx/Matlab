function solvetype = z_setsolve(zchan,n_surf,code,pnum,varargin)
% retval = z_setsolve(zchan,n_surf,code,pnum,solvetype,param #1, param #2, ...)
%
% This item returns data about any solve on any surface. The syntax is 
% GetSolve,surface,code
% 
% where code is an integer code indicating which surface parameter the solve
% data is for. The solve data is returned in the following
% formats, depending upon the code value.
% 
% if code is parameter or extra data value, pnum indicates which parameter
% or extra data value, otherwise, pnum is not used.
% 
% The solvetype is an integer code, and the parameters have meanings that
% depend upon the solve type; see
% the chapter Solves for details. See also SetSolve.
% 
% GetSolve Code Returned data format
% 0, curvature solvetype, parameter1, parameter2
% 1, thickness solvetype, parameter1, parameter2, parameter3
% 2, glass solvetype, pickupsurf
% 3, semi-diameter solvetype, pickupsurf
% 4, conic solvetype, pickupsurf
% 5-12, parameters 1-8 solvetype, pickupsurf, offset, scalefactor
% 1001+, extra data values
% 1+
% solvetype, pickupsurf, scalefactor

% validate code
if ischar(code)
    switch lower(code)
        case {'curv','curvature'}
            code = 0;
        case {'thick','thickness'}
            code = 1;
        case {'glass'}
            code = 2;
        case {'semi-diameter','sdia'}
            code = 3;
        case {'conic'}
            code = 4;
        case {'parameter','parm','param'},
            code = 4 + pnum;
        case {'extra data','extra'},
            code = 1000 + pnum;
        otherwise
            error(['code ' code ' is not recognized code type']);
    end
else
    if code < 0,
        error([num2str(code) ' is not a valid code']);
    end
end

cmdstr = ['SetSolve,' num2str(n_surf) ',' num2str(code) ];
for ii = 1:length(varargin),
    cmdstr = [cmdstr ',' num2str(varargin{ii})];
end

solvetype = ddereq(zchan,cmdstr,[1 1]);
