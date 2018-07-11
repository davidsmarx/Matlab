% rays = mexReadZRD(buffer, sOptions);
%
% buffer is an unit8 array of bytes read from the binary ZRD file
% sOptions = struct(...
%   'SourceObject', object number of source from zemax model
%   ,'DectectorObject', object number of detector from zemax model
%   );
% 
% rays = array of struct:
% % 	const char *field_names[] = { "status"
% 	,"level"
% 	,"hit_object"
% 	,"hit_face"
% 	,"unused"
% 	,"in_object"
% 	,"parent"
% 	,"storage"
% 	,"xybin", "lmbin"
% 	,"index", "starting_phase"
% 	,"x", "y", "z"
% 	,"l", "m", "n"
% 	,"nx", "ny", "nz"
% 	,"path_to", "intensity"
% 	,"phase_of", "phase_at"
% 	,"exr", "exi", "eyr", "eyi", "ezr", "ezi"
%   ,"wavelength"};
