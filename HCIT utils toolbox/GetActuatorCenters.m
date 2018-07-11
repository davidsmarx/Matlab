function [XAused, YAused, Sdmind] = GetActuatorCenters(base_static)
    % base_static is from python model config, load it from saved matfile
    %
    % (XAused, YAused) are list of actuator centers in same space as image
    % pixels
    
    
    if nargin == 0,
        base_static = struct(...
            'infmatfn','/home/bseo/HCITtest/trunk/HLC_alignment/Support/infmat.hlc.ideal.48.1.fits' ... % (i,) for i in [1,2]],
            ,'ijlistfn','/home/bseo/HCITtest/trunk/HLC_alignment/Support/ijlist_hlc_4608_all_AFTA.fits' ...
            ,'xycornfn','/home/bseo/HCITtest/trunk/HLC_alignment/Support/xycorn.hlc.ideal.48.%d.fits'... % (i,) for i in [1,2]],
            );
    end
                                   
    %%% ijlist is a list of actuators that are used
    if isfield(base_static,'ijlist'),
        ijlist = base_static.ijlist;
    elseif isfield(base_static,'ijlistfn'),
        ijlist = fitsread(PathTranslator(base_static.ijlistfn));
    else
        error('no ijlist');
    end
    
    % python: i,j = ij % self.nact[idm], ij//self.nact[idm]  % 0-offset
    % (should get # actauators = 48 from fitsinfo)
    NxAct = 48; NyAct = 48;
    i = mod(ijlist,NxAct); j = floor(ijlist/NyAct);
    [idm1, jdm1] = filterdata(j<NyAct, i+1, j+1);
    [idm2, jdm2] = filterdata(j>=NyAct, i+1, j-NxAct+1);
    % figure, plot(idm1, jdm1, 'o'), grid, title('DM1 i j list'), xlabel('i'), ylabel('j')
    Ndm1 = length(idm1); % # of dm1 actuators in the ijlist

    % output ijlist for use with output Act_ht
    Sdmind = struct(...
        'NxAct', NxAct ...
        ,'NyAct', NyAct ...
        ,'ijlist', ijlist ...
        ,'idm1', idm1, 'jdm1', jdm1 ...
        ,'idm2', idm2, 'jdm2', jdm2 ...
        ,'Nuse', Ndm1 ...
        );

    % actuator grid (should get # actauators = 48 from fitsinfo)
    xa = (-24:23)'*7; % 7 pixels per actuator
    % coordinates of actuators in the ijlist
    [XAused, YAused] = deal(zeros(size(idm1)));
    for ij = 1:Ndm1,
        XAused(ij) = xa(idm1(ij));
        YAused(ij) = xa(jdm1(ij));
    end
    % (XAused, YAused) are list of actuator centers in same space as image
    % pixels

end % GetActuatorCenters
