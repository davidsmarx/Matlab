function winpath = PathTranslator(s383path)
% winpath = PathTranslator(s383path)

if ~ispc, %isequal(computer, 'GLNXA64'),
    winpath = s383path;
    return
end

% init output
winpath = '';

A = regexp(s383path,'/','split');

if length(A) <= 2,
    winpath = s383path;
    return
end

if strcmp(A{1}, '~'),
    winpath = ['Y:\' A{2} '\' A{3:end}];
    
elseif ~strcmp(s383path(1), '/')
    winpath = s383path; % it's a relative path
    
elseif strcmp([A{2} '/' A{3}], 'home/bseo'),
    winpath = strjoin({'Y:\bseo', A{4:end}}, '\');
    
elseif strcmp([A{2} '/' A{3}], 'home/dmarx'),
    winpath = strjoin({'Y:', A{4:end}}, '\');
    
elseif isequal([A{2} '/' A{3}],'proj/hcit2')
    % \\s383-nfs = /proj
    %winpath = '\\s383-nfs\hcit2'; %
    winpath = ['Y:\' strjoin({A{5:end}}, '\')]; % A{4} = dmarx
    
elseif isequal([A{2} '/' A{3} '/' A{4} '/' A{5}],'proj/hcit/home/dmarx')
    % \\s383-nfs = /proj
    %winpath = '\\s383-nfs\hcit'; %
    winpath = ['Y:\' strjoin({A{6:end}}, '\')];
    
elseif isequal([A{2} '/' A{3}],'proj/afta-im')
    % \\s383-nfs = /proj
    winpath = '\\s383-nfs\afta-im'; %

elseif isequal([A{2} '/' A{3} '/' A{4}],'proj/hcit/data')
    winpath = strjoin({'Y:\ln_hcit_data', A{5:end}}, '\');

elseif isequal([A{2} '/' A{3} '/' A{4}],'proj/mcb/data')
    % \\s383-nfs = /proj
    %winpath = strjoin({'\\muscle5.jpl.nasa.gov\mcb', A{4:end}}, '\');
    winpath = strjoin({'Y:\ln_mcb_data', A{5:end}}, '\');

elseif isequal([A{2} '/' A{3}],'net/spud-data')
    % mcb camera data
    winpath = 'Y:\HCIT\ln_spud_data_Data';
    A = {A{1:3},A{5:end}}; % remove /Data/

elseif isequal([A{2} '/' A{3}],'net/piaa-data') || isequal([A{2} '/' A{3}],'proj/piaa-data')
    % mcb camera data
    winpath = strjoin({'X:', A{5:end}}, '\');
    %A = {A{1:3},A{5:end}}; % remove /Data/

elseif isequal([A{2} '/' A{3}], 'proj/piaacmc')
    % piaacmc share and camera images
    winpath = strjoin({'W:', A{4:end}}, '\');
    
elseif length(A) >= 4 && isequal([A{2} '/' A{3} '/' A{4}],'proj/dst/data')
    winpath = strjoin({'Y:\ln_dst_data', A{5:end}}, '\');
    
elseif isequal([A{2} '/' A{3}], 'proj/piaacmc'),
    % piaacmc
    winpath = strjoin({'W:', A{4:end}}, '\');
    
elseif isequal(A{1}, '..') || isequal(A{1}, '.'),
    % relative path, no need to translate
    winpath = s383path;
    
else,
    warning(['no translation for ' s383path]);
    winpath = s383path;
    return
end

% for ii = 4:length(A),
%     winpath = [winpath '\' A{ii}];
% end
