function [SMaps] = ZRD_ReadZRD(filename, sOptions)
% [rays] = ZRD_ReadZRD(filename, sOptions)
%
% sOptions = struct( ...
%   'SourceObject', 54 ...
%  ,'DectectorObject', 46 ...
%   );

if nargin == 0 || isempty(filename),
    [fname, pathname] = uigetfile({'*.zrd'},'Select ZRD file');
    if isequal(fname,0), return, end
    filename = [pathname fname];
end

% get the size of the file in bytes
S = dir(filename);

% open the binary file and read the whole thing into a buffer of bytes
fid = fopen(filename,'r');
buffer = fread(fid, S.bytes, '*uint8');

%
fclose(fid);

% check the size of the buffer
buffersize = length(buffer);
if buffersize ~= S.bytes, error(['read only ' num2str(buffersize) ' bytes']); end

% sort out the rays in a mex routine
SMaps = mexReadZRD(buffer, sOptions);
for isou = 1:length(SMaps),
    SMaps(isou).Phase = angle(SMaps(isou).Field);
    %CohIrr = IrrMap .* abs(FldMap).^2 ./ AmpMap.^2;
    SMaps(isou).CohIrr = SMaps(isou).Field.*conj(SMaps(isou).Field);
end

end % main
