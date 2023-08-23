function writeyaml(fn, data)
% writeyaml(fn, data)
%
% my overlay of this yaml module to write yaml files the way I like

fid = fopen(fn, 'wt');
fprintf(fid, yaml.dump(data, "block"));
fclose(fid);

