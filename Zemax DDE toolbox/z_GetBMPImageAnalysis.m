function S = z_GetBMPImageAnalysis(zchan, textfilename)
% S = z_GetBMPImageAnalysis(zchan, textfilename)
% 
% calls ZEMAX Bitmap image analysis
% the settings are the default settings for the lens
%
% Results are stored in ASCII form in textfilename, the
% default is [pwd '\imageanalysis.txt']

% efficiency = fraction of rays to reach the image plane
% img = [x,y] location of rays in the image plane
% field = [x,y] location of rays at source
% pupil = [x,y] location of rays at the pupil
% wave = wavelength# for each ray

% validate inputs
textfilename = [pwd '\bmpimageanalysis.txt'];
% file for zemax settings
settingsfilename = [pwd '\IbmSettings.cfg'];
% flag determines use/save settings. see p. 634 in Zemax manual
settingsflag = '0'; % 0 => use default settings, and save the settings used to settingsfilename

% Geometric Image Analysis
typecode = 'Ibm';

% make the call
cmdstr =...
    ['GetTextFile, ' textfilename ', ' typecode ', ' settingsfilename ', ' settingsflag];
retval = ddereq(zchan,cmdstr,[1 1]);

% wait until zemax is ready
while ~ddereq(zchan,'GetVersion',[1 1]), pause(0.25); end

% get the image data from the textfile
%S = ReadImageanalysisTxt(textfilename);

% for now just read the image file written in C:\Program Files\ZEMAX\IMAFiles
S = imread(['C:\Program Files\ZEMAX\IMAFiles\BitmapImage.JPG']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = ReadImageanalysisTxt(textfile)

global UM MM;

fid = fopen(textfile,'rt');
if fid == -1, error(['error opening text file ' textfile]); end

% read in all of the ray data
C = textscan(fid,'%d %f %f %f %f %f %f %f %f','headerlines',12);
[rayno, xfield, yfield, px, py, wave, wgt, xi, yi] = deal(C{:});
field = [xfield yfield]*MM;
pupil = [px py]*MM;
img   = [xi yi]*MM;

% skip to last line to get efficiency value
if fseek(fid,-52,'eof') == -1, error(ferror(fid)); end
efficiency = fscanf(fid,'%*s %f');
efficiency = 0.01*efficiency; % convert percent to fraction

fclose(fid);

end
