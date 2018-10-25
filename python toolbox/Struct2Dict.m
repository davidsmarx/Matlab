function adict = Struct2Dict(astruct)
% 

% matrices with ndims >= 2 must be handles separately
fnames = fieldnames(astruct);
for ii = 1:length(fnames),
    switch class(astruct.(fnames{ii})),
        case {'double','int8','int16','int32','int64'},
            atmp = astruct.(fnames{ii});
            astruct.(fnames{ii}) = Mat2npArray(atmp);
        case 'struct',
            atmp = Struct2Dict(astruct.(fnames{ii}));
            astruct.(fnames{ii}) = atmp;
        case 'logical'
            atmp = Mat2npArray(astruct.(fnames{ii}));
            astruct.(fnames{ii}) = atmp;
            
        otherwise,
            % do nothing
    end % switch class
end

adict = py.dict(astruct);
