function simparms = SIM_wideangle_grid(sp)
% simparms = SIM_wideangle_grid(sp)
%
% sets nfov, fovx, fovy, fov fields in the simparms struct to correspond to
% the 57 pre-defined field points.
    
global P;

simparms = sp;

simparms.nfov = 57;

xtmp = [1:2:7 1:2:7 1 3 1 3 4+1/3 5+2/3]*P;
ytmp = [[1 1 1 1] [3 3 3 3] [5 5] [7 7] 5+2/3 4+1/3]*P;
simparms.fovx = [xtmp  xtmp -xtmp -xtmp 0]';
simparms.fovy = [ytmp -ytmp  ytmp -ytmp 0]';
simparms.fov  = max(sqrt(simparms.fovx(:).^2 + simparms.fovy(:).^2));

