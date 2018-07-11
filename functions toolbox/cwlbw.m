function [cwl, bw] = cwlbw(lam,td,amp,varargin)
% [cwl, bw] = cwlbw(lam,td,amp,options)
% 
% inputs:
%   lam [um]
%   td  amplitude [dB] or linear (default), if td is a matrix,
%       then number of rows of td must equal length(lam),
%       and cwl and bw is found for each column
%   amp threshold for calculating BW [same units as td]
% output:
%   cwl is mean [um]
%   bw is full bandwidth [um]
% options:
%   'dB' if td is in [dB]

if nargin == 0, error('usage: [cwl, bw] = cwlbw(lam,td,amp,options)'); end

if isempty(strmatch('dB',varargin)), td = abs(td); end

% check whether td is vector or matrix
[nr, nc] = size(td);
if any([nr, nc] == length(td(:))), % td is a vector
   td = td(:);
   nc = 1;
else, % td is a matrix
   if nr ~= length(lam), error('nrows of td must = length(lam)'); end
end

for iic = 1:nc,
   tdi = td(:,iic);
   
   ii = [1:length(lam)];
   na = tdi > amp;
      
   ia = ii(na);
   if isempty(ia),
      if abs(tdi(1))>abs(tdi(end)),
         cwl(iic) = -1000;
      else,
         cwl(iic) = 1000;
      end
      bw(iic) = 1000;
      break;
   end
   
   i1 = ia(1); i2 = ia(end);
   
   if i1 == 1,
      cwl(iic) = -1000;
      bw(iic)  = 1000;
      break
   elseif i2 == ii(end),
      cwl(iic) = 1000;
      bw(iic) = 1000;
      break
   end
   
   a1 = interp1(tdi(i1-1:i1),lam(i1-1:i1),amp,'spline');
   a2 = interp1(tdi(i2:i2+1),lam(i2:i2+1),amp,'spline');
   
   cwl(iic) = mean([a1 a2]);
   bw(iic)  = a2-a1;
   
end % for each column

return