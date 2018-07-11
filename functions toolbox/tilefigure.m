% TILEFIGURE
%     This function allows you to take N saved figures and tile them 
%     into one figure.  The inputs are the filenames of the saved
%     figures (.m extension not required).  The last input is the 
%     filename of the tiled plot.  
%
%     EXAMPLE
%       >> tilefigure('foo.m','bar','foobar')
%
%     Figures with text other than xlabels, ylabels, zlabels and 
%     short titles are not recommended as this tends to crowd the 
%     tiled figures.  Figures with subplots are not compatible.  
%
%     Tiling begins in the bottom left hand corne and moves right then up. 
% 
%
%     You can add to a tile in edit mode:
%     EXAMPLE 
%       >> Do you wish to add to figure <file>? (y/n): y
%       :colorbar             (adds a colorbar to side of figure)
%       :text(0,0,'Origin')   (adds text to (0,0) of <file>)
%       :-1                   (quits edit mode)
%       >>
%     Some exceptions are plot, plot3, subplot and figure filenames.
%     You will be making the additions blind so keep the other figures 
%     up.  To quit edit mode type -1.
%
%     Requirements:  A .m and .mat file must exists for each figure.
%
%     Note:  Remember to surrond each input with single quotes (').
%     
%     Andrew Wheatley
%     awheat@irus.rri.on.ca
%     Tested on Matlab 5.0.0.4064

  function tilefigure(varargin)

% check number of input variables
 
  n = nargin;
  if n < 3
    str = ['A minimum of 2 saved files plus a new filename are required.'];
    disp(str)
    return
  end

% test for overwrite
% by opening last input with read only permission

% convert cell varargin to readable struct field 
  field = ['filename'];
  tilefile = cell2struct(varargin(n),field,2);
  tilefile = tilefile.filename;

% check for .m extension
  if length(tilefile) > 2 & ...
   strncmp(tilefile(length(tilefile)-1:length(tilefile)),'.m',2) != 1
    tilefile(length(tilefile)+1:length(tilefile)+2) = '.m';
  end
 
% open the new file for reading 
  newfid = fopen(tilefile,'r');
  if newfid != -1
    % if file exists query user on overwrite permission
      str = ['Overwrite ',tilefile, '? [Y/N]: '];
      i = input(str,'s');
    
    % check user input
      if isempty(i) | i == 'n' | i == 'N'
        str = ['New filename: '];
        j = input(str,'s');
        varargin(n) = cellstr(j);
      
        % convert cell varargin to readable struct field 
          field = ['filename'];
          tilefile = cell2struct(varargin(n),field,2);
          tilefile = tilefile.filename;

        % check for .m extension
          if length(tilefile) > 2 & ...
           strncmp(tilefile(length(tilefile)-1:length(tilefile)),'.m',2) != 1
            tilefile(length(tilefile)+1:length(tilefile)+2) = '.m';
          end

          fclose(newfid);
      else      
        fclose(newfid);
      end
  end
  
% make the dimension for the subplot as square as possible
% nrow is number of rows
% ncol is number of columns
% preference is given to rows

  nrow = floor(sqrt(n-1));

  if nrow*nrow == n-1          % if n-1 is a square
    ncol = nrow;
  elseif nrow*(nrow+1) >= n-1  % else n(n-1)
    ncol = nrow;
    nrow = nrow+1;
  else                         % otherwise n^2 is best approximation
    nrow = nrow+1;
    ncol = nrow;
  end

% initialize the left,bottom coordinates of the first tile
  initial_left = 0.12;
  initial_bottom = 0.12;

% define the usable space in the new figure
  usable_width = 0.9;
  usable_height = 0.9;

% define percentage of tile that will make up each figure 
% remaining percentage of tile is used as space between tiles
  fig_width = 0.75;
  fig_height = 0.70;

% calculate the tile width and height
  sub_width = usable_width*(fig_width/ncol);
  sub_height = usable_height*(fig_height/nrow);

% initialize row and column counters used for tiling
  rowcnt = 1;
  colcnt = 1;

% open the new file for reading and writing 
  newfid = fopen(tilefile,'w+');
  if newfid == -1
      str = ['Cannot open ', tilefile, '.'];
      disp(str)
      return
  end

% load in the first file
% convert cell varargin to readable struct field 
  field = ['filename'];
  file1 = cell2struct(varargin(1),field,2);
  file1 = file1.filename;
  
% check for .m extension
  if length(file1) > 2 & ... 
   strncmp(file1(length(file1)-1:length(file1)),'.m',2) == 1
    file1 = file1(1:length(file1)-2);
  end

% write the comments and initial handle that will initiate the new file
  fprintf(newfid, '%% %s\n', tilefile);
  fprintf(newfid, '%% Tiled figure containing the plots');
  for i = 1:n-1
    field = ['filename'];
    file = cell2struct(varargin(i),field,2);
    if i < n-1
      fprintf(newfid, ' %s,', file.filename);
    else
      fprintf(newfid, ' and %s.', file.filename);
    end
  end

  fprintf(newfid, '\n\nload %s\n\n', file1);
  fprintf(newfid, 'a = figure(\''Color\'',[0.8 0.8 0.8], ...\n');
  fprintf(newfid, '           \''Colormap\'',mat0, ...\n');
  fprintf(newfid, '           \''Position\'',[254 253 512 384]);\n');
 
% now loop through the first n-1 inputs 
% and write to the new file (input n)
  
  for i = 1:n-1
    % calculate left, bottom coordinates of tile i
      left = initial_left + usable_width*(colcnt-1)/ncol;
      bottom = initial_bottom + usable_height*(rowcnt-1)/nrow;

    % convert cell to struct
      field = ['filename'];
      file = cell2struct(varargin(i),field,2);
      file = file.filename;

    % check for .m extension
      if length(file) > 2 & ...
       strncmp(file(length(file)-1:length(file)),'.m',2) != 1
        file(length(file)+1:length(file)+2) = '.m';
      end

    % open the each file with read only permission 
      oldfid = fopen(file,'r');
      if oldfid == -1
        str = ['Cannot find ', file, '.'];
        disp(str)
        return
      end
  
    % index through each line in the file until you find the first one
    % that begins with 'b = axes'
      index = [];
      while feof(oldfid) != 1
        buffer = fgets(oldfid);
        index = find(buffer == '=');

        % case of more than one instance of '=' on a line
          index = min(index); 
          if strncmp(buffer(index-2:index+5),'b = axes',8) == 1
            break;
          end
      end

      fprintf(newfid, '%s', buffer); % 'b = axes' line

    % begin writing oldfid to newfid
      while feof(oldfid) != 1
        buffer = fgets(oldfid);
         
        % check for another instance of b = axes in old file
        % if there is then we have a subplot and must quit
          if strncmp(buffer(1:8),'b = axes',8) == 1
            str = ['Subplot figure error in ', file,'.'];
            disp(str)
            return
          end

          if length(buffer) > 9 & strncmp(buffer(3:12),'ColorOrder',10)
            % write the position of the tile to newfid after 'ColorOrder'
              fprintf(newfid,'    \''Position\'', [%g %g %g %g], ...\n', ...
                      left, bottom, sub_width, sub_height);
          else
            % otherwise just keep on writing
              fprintf(newfid, '%s', buffer);
          end
      end  % end while

    % close oldfid
      fclose(oldfid);

    % give user oppurtunity to add to tile i
      str = ['Do you wish to add anything to ', file, '? [Y/N]: '];
      j = input(str,'s');
      if ~isempty(j) & (j == 'y' | j == 'Y')
        disp('Entering edit mode: (-1 to quit)')
        edit_mode = 1;
        while edit_mode == 1
          j = input(':','s');
          if strncmp(j(1:2),'-1',2) == 1
            % quits edit mode
              edit_mode = 0;
          else
            fprintf(newfid,'\n%s\n',j);
            edit_mode = 1;
          end % end inner if

        end % end while

      end % end outer if
      
    % each element can be thought of as a matrix element
    % the following 8 lines determines the element label
    % ie element(rowcnt,colcnt)
      colcnt = colcnt+1;
      if colcnt > ncol
          colcnt = 1;
          rowcnt = rowcnt+1;
          if rowcnt > nrow
              rowcnt = 1;
          end
      end

    % load in the next file
      if i < n-1
          field = ['filename'];
          file = cell2struct(varargin(i+1),field,2);
          file = file.filename;
          % check for .m extension
            if strncmp(file(length(file)-1:length(file)),'.m',2) == 1
              file = file(1:length(file)-2);
            end

          % begin next file
            fprintf(newfid, '\nload %s\n\n', file);
      end
  end % end for loop

% close the new file
  fclose(newfid);

% tell user what to do
% check for .m extension
  if strncmp(tilefile(length(tilefile)-1:length(tilefile)),'.m',2) == 1
    tilefile = tilefile(1:length(tilefile)-2);
  end

  str = ['To see new figure, type "', tilefile,'".'];
  disp(str)

