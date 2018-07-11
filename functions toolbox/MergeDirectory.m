function emsg = MergeDirectory(dirFrom, dirTo)
% 

emsg = '';

if ~exist('dirFrom','var') && ~isempty('dirFrom'),
    dirFrom = uigetdir('C:\Users\dmarx\inSync Share','Source Directory');
end
if ~exist('dirTo','var') && ~isempty('dirTo')
    dirTo = uigetdir('C:\users\dmarx','Destination Directory');
end

sDirFrom = dir(dirFrom);
sDirTo   = dir(dirTo);

for ii = 1:length(sDirFrom),
    
    rval = 0;
    
    if any(strcmp(sDirFrom(ii).name,{'.','..'})),
        continue;
    end
    
    if sDirFrom(ii).isdir,
        continue;
    end
    
    jj = find(strcmp(sDirFrom(ii).name, {sDirTo.name}));
    
    if isempty(jj)
        rval = dos(['copy "' dirFrom '\' sDirFrom(ii).name '" "' dirTo '"'])
        

    elseif sDirFrom(ii).datenum > sDirTo(jj).datenum,
        % ~isempty, choose more recent
        rval = dos(['copy ' dirFrom '\' sDirFrom(ii).name ' ' dirTo])
    end
    
    if rval ~= 0,
        disp('failed:');
        disp(['copy "' dirFrom '\' sDirFrom(ii).name '" "' dirTo '"']);
    end
    
end
