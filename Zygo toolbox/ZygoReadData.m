function [phase, header, intens, success] = ZygoReadData(filename)
%
% [phase] = ZygoReadData(filename)
% [phase, header, intens, success] = ZygoReadData(filename) 
%
% Read all the data from a Zygo "dat" file format.
%
% Input:
%   filename = the name of the Zygo file to read.
%
% Output:
%   phase   = phase data array in the same units as header.Wavelengthin
%             (usually meters)
%   header  = a structure containing selected contents from the file header.
%   intens  = an array containing the intensity data from the file.
%   success = a scalar indicating success (1) or failure (0) of the function.
%
% Versions:
%   09/14/98 Ladd Wheeler  Original version
%   08/30/01 Cathy Ohara   Scale phase data and also return intensity.
%   01/10/06 Cathy Ohara   Use coordinates in header (PhaseOriginX,
%                          PhaseOriginY, PhaseWidth, PhaseHeight) to register
%                          the phase map.  Also pass out "phase" first.
%   07/02/15 D. Marx       Copied to my toolbox, added uigetfile

phase   = [];
header  = [];
intens  = [];
success = 0;

nout = nargout;

if nargin == 0,
    [fname, pathname, fi] = uigetfile({'*.dat'},'Select Zygo Data File');
    filename = [pathname fname];
end

% Open the file for read-only binary access with big-endian format.
file = -1;
[file, message] = fopen (filename, 'r', 'ieee-be');
if file == -1
    error (message);
end

% Now read the header data.
[header, success, offset] = readzygoheaderdata (file);
if success == 0
    error ('Failure in reading Zygo header data');
end

% Now read the phase data.
if (nout > 2)
   [phase, intens, success] = readzygophasedata (file, header, offset);
else
   [phase] = readzygophasedata (file, header, offset);
   success = 1;   %CAT:  cheesy but speedy
end

if success == 0
    error ('Failure in reading Zygo phase data');
end

% Close the file.
fclose (file);


% Use the header information (PhaseOriginX,PhaseOriginY,PhaseWidth,PhaseHeight)
% to reconstruct the phase map.

%  oo = zeros(1024);

  x = header.PhaseOriginX;
  y = header.PhaseOriginY;
  w = header.PhaseWidth;
  h = header.PhaseHeight;
  imw = header.IntensWidth;
  imh = header.IntensHeight;

  oo = zeros(imh, imw);
  
  %oo(y+1:(y+h),x+1:(x+w)) = phase;
  oo = phase;
  oo = oo * header.Wavelengthin;
  
  phase  = oo;
  status = 1;
  
return

% --------------------------------------------------------------------------

function [header, success, header_size] = readzygoheaderdata (file)
%
% ReadZygoHeaderData reads selected data from the specified Zygo file.
%
% Input:
%   file = the file id of the Zygo file to read.
%
% Output:
%   header = a structure containing selected contents from the file header.
%   success = a scalar indicating success (1) or failure (0) of the function.
%
% Versions:
%   09/15/98 Ladd Wheeler  Original version

success = 1;

% Read Magic Number
status = fseek (file, 0, 'bof');
if  status ~= 0
    success = 0;
    error ('Failure to seek to 0 in file.');
end
[magic_number] = fread (file, 1, 'int32');
%sprintf(int2str(magic_number))
%if  magic_number ~= -2011495569
%    success = 0;
%    error ('Magic number in file is not valid.');
%end

% Read Header Format
header_format = fread (file, 1, 'int16');
%sprintf(int2str(header_format))
%if  header_format ~= 1
%    success = 0;
%    error ('Header format value in file is not 1.');
%end

% Read Header Size
header_size = fread (file, 1, 'int32');
%sprintf(int2str(header_size))
%if  header_size ~= 834
%    success = 0;
%    error ('Header size value in file is not 834.');
%end

% Read Intensity Description Data
status = fseek (file, 48, 'bof');
if  status ~= 0
    success = 0;
    error ('Failure to seek to Intensity Origin values in file.');
end
header.IntensOriginX = fread (file, 1, 'int16');
header.IntensOriginY = fread (file, 1, 'int16');

header.IntensWidth = fread (file, 1, 'int16');
header.IntensHeight = fread (file, 1, 'int16');

header.NBuckets = fread (file, 1, 'int16');
header.IntensRange = fread (file, 1, 'uint16');
header.IntensBytes = fread (file, 1, 'int32');

% Read Phase Description Data
header.PhaseOriginX = fread (file, 1, 'int16');
header.PhaseOriginY = fread (file, 1, 'int16');

header.PhaseWidth = fread (file, 1, 'int16');    
header.PhaseHeight = fread (file, 1, 'int16');    %Byte #70 & #71

header.PhaseBytes = fread (file, 1, 'int32');

% Read Time Stamp and Comment
header.TimeStamp = fread (file, 1, 'int32');
header.Comment = fread (file, 82, 'char');
header.Comment = char(header.Comment)';

% Read Other
header.Source = fread (file, 1, 'int16');
% ***************************************************####################
%dummy = fread (file, 1, 'char');
% ***************************************************####################
header.IntfScaleFactor = fread (file, 1, 'float32');    %Scale factor, Byte 164+
header.Wavelengthin = fread (file, 1, 'float32');
header.NumericAperature = fread (file, 1, 'float32');
header.ObliquityFactor = fread (file, 1, 'float32');    %Obliquity, Byte 176+
skip = fread (file, 1, 'float32');
header.CameraRes = fread (file, 1, 'float32');

% Read Next Clumps of Descriptives
status = fseek (file, 218, 'bof');
if  status ~= 0
    success = 0;
    error ('Failure to seek to PhaseRes in file.');
end
header.PhaseRes = fread (file, 1, 'int16');            %Phase Resolution, Byte 218+
header.MinimumAreaSize = fread (file, 1, 'int32');
header.DisconAction = fread (file, 1, 'int16');
header.DisconFilter = fread (file, 1, 'float32');
header.ConnectionOrder = fread (file, 1, 'int16');
header.DataSign = fread (file, 1, 'int16');
header.CameraWidth = fread (file, 1, 'int16');
header.CameraHeight = fread (file, 1, 'int16');

status = fseek (file, 300, 'bof');
if  status ~= 0
    success = 0;
    error ('Failure to seek to PhaseAvgs in file.');
end
header.PhaseAvgs = fread (file, 1, 'int16');
header.SubtractSysErr = fread (file, 1, 'int16');

status = fseek (file, 364, 'bof');
if  status ~= 0
    success = 0;
    error ('Failure to seek to RemoveTiltBias in file.');
end
header.RemoveTiltBias = fread (file, 1, 'int16');
header.RemoveFringes = fread (file, 1, 'int16');
header.MaxAreaSize = fread (file, 1, 'int32');
skip = fread (file, 1, 'int16');
header.PreConnectFilter = fread (file, 1, 'float32');

% --------------------------------------------------------------------------

function [phase, intens, success] = readzygophasedata (file, header, offset)
%
% ReadZygoPhaseData reads Phase data from the specified Zygo file.
%
% Input:
%   file = the file id of the Zygo file to read.
%   header = a structure containing selected contents from the file header.
%
% Output:
%   phase = an array containing the phase data.
%   success = a scalar indicating success (1) or failure (0) of the function.
%
% Versions:
%   09/15/98 Ladd Wheeler  Original version.

success = 1;

nout = nargout;

% Check if any Phase Data is in file
if  header.PhaseBytes == 0
    error ('There is no Phase data in the file.');
end

% Try to read intensity data
if (nout > 1)
%offset = 834;
status = fseek (file, offset, 'bof');
if  status ~= 0
    success = 0;
    error ('Failure to seek the Intensity Data in the file.');
end
intens = fread (file, [header.IntensWidth, header.IntensHeight], 'int16')';
end

%   intens = zeros (header.IntensHeight, header.IntensWidth);
%   for row = 1:header.IntensHeight
%       for col = 1:header.IntensWidth
%           intens(row,col) = fread (file, 1, 'int16');
%       end
%   end

% Seek to the beginning of the Intensity Data
%offset = 834 + header.IntensBytes;
offset = offset + header.IntensBytes;
status = fseek (file, offset, 'bof');
if  status ~= 0
    success = 0;
    error ('Failure to seek to the Phase Data in the file.');
end

% Preallocate the phase array to speed up the execution.

% Read the Phase Data
%   According to Zygo documentation, the data is written in "row-major order".
%   This unit assumes that the scan starts at the upper left (row 1, col 1)
%   and moves left-to-right along EACH row.

phase = fread (file, [header.PhaseWidth, header.PhaseHeight],'int32')';

%phase = zeros (header.PhaseHeight,  header.PhaseWidth);
%for row = 1:header.PhaseHeight
%    for col = 1:header.PhaseWidth
%        phase(row,col) = fread (file, 1, 'int32');
%    end
%end

% Multiply by a bunch of weird scale factors.  See Andrew's IDL routine, zygo.pro.
datasign   = header.DataSign; datasign=1-2*datasign;
resolution = 4096 * 8^(header.PhaseRes);
mask       = phase < 2147483640;

phase = (datasign*header.IntfScaleFactor*header.ObliquityFactor/resolution) * phase .* mask;
%phase = (datasign*header.ObliquityFactor/resolution) * phase .* mask;
%disp('CAT:  Not setting background to minimum in readzygophasedata.m') 
%phase = phase + min(phase(:))*(mask==0);
