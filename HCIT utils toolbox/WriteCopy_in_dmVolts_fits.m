function [dm1, dm2] = WriteCopy_in_dmVolts_fits(dm1, dm2,...
    dmfitssrc, dmfitsdest, listKeyVals)
% [dm1, dm2] = WriteCopy_in_dmVolts_fits(dm1, dm2, dmfitssrc, dmfitsdest, listKeyVals)
%
% copy all the keywords for dmfitssrc
% replace dm volts with inputs
% dm1, dm2:
%     empty => use dm volts from dmfitssrc
%     48x48 array of Volts => write arrays directly to new dmfitsdest
%     fits_fn => read dm volts from fits_fn
%     struct('dm_displacement', dm ...
%           ,'gain_fn','/proj/mcb/data/dB_dm_gain/dm481.gain.fits')
%           => use gain to convert displacement to volts, and write volts
%           to dmfitsdest
%
% if listKeyVals is a cell array of keyword, value pairs, update the
% keywords with the new values.
%
% examples:
%
% WriteCopy_in_dmVolts_fits(30*ones(48), 30*ones(48),...
%     PathTranslator('/proj/mcb/data/EFC/HLC_model/run470/in/it00000dm.fits'),...
%     'results\run000\in\it00000dm.fits');
%
% [dm1, dm2] = WriteCopy_in_dmVolts_fits(...
%      PathTranslator('/proj/mcb/data/dB_PR_Kern/gsomc_no419_dm1.fits'),...
%      PathTranslator('/proj/mcb/data/dB_PR_Kern/gsomc_no419_dm2.fits'),...
%      PathTranslator('/proj/mcb/data/EFC/HLC_model/run470/in/it00000dm.fits'),...
%      'results\run000\in\it00000dm.fits',{'ICOLOR0',2});

import matlab.io.*

finfo = fitsinfo(dmfitssrc);

dm1 = ValidateDMinput(dm1, dmfitssrc, 1);
dm2 = ValidateDMinput(dm2, dmfitssrc, 2);

% create new fits file
% prefix '!' means overwrite if already exists
fptr = fits.createFile(['!' dmfitsdest]);

% update keyword values
if exist('listKeyVals','var') && ~isempty(listKeyVals) && iscell(listKeyVals),
    [nkeys, n2] = size(listKeyVals);
    for ik = 1:nkeys,
        finfo.PrimaryData.Keywords = FitsSetKeywordVal(finfo.PrimaryData.Keywords, ...
            listKeyVals{ik,1}, listKeyVals{ik,2});
    end
end

% if dm1, dm2 inputs are displacements, use dm gain to convert to Volts


% primary HDU
writeHDU(fptr, dm1, finfo.PrimaryData.Keywords);

% Image HDU
writeHDU(fptr, dm2, finfo.Image.Keywords);

fits.closeFile(fptr);

end % main

function writeHDU(fptr, dm, Keywords)
        
import matlab.io.*

    fits.createImg(fptr,'double',size(dm));
    fits.writeImg(fptr,dm);
    
    [~, ikExt] = FitsGetKeywordVal(Keywords, 'EXTEND');
    [~, ikEnd] = FitsGetKeywordVal(Keywords, 'END');
    
    if isempty(ikExt) || isempty(ikEnd),
        return
    end
    
    ChangeToInt('NCOLOR');
    ChangeToInt('NPAIRS');
    ChangeToInt('XSIZE');
    ChangeToInt('YSIZE');
    ChangeToInt('NEXP');
    ChangeToInt('NUMIM');
    ChangeToInt('ICOLOR0');
    ChangeToInt('ICOLOR1');
    ChangeToInt('ICOLOR2');
    ChangeToInt('ICOLOR3');
    ChangeToInt('ICOLOR4');
    ChangeToInt('PROBEF');
    
    
    function ChangeToInt(key)
        [val, ik] = FitsGetKeywordVal(Keywords, key);
        if ~isempty(val), Keywords{ik,2} = int32(val); end
    end
    
    for ikey = (ikExt+1):(ikEnd-1),
        ktmp = Keywords(ikey,:);
        %if any(cellfun(@isempty,ktmp{1:2})),
        if isempty(ktmp{2})
            ktmp{2} = 0;
            disp('Empty Keyword entry:');
            disp(ktmp);
        end
        %disp(ktmp);
        fits.writeKey(fptr,ktmp{1:2});
    end
        
end % writeHDU

function dm = ValidateDMinput(dm_input, dmfitssrc, dm1ordm2)

    U = CConstants;
    
    if isempty(dm_input),
        if dm1ordm2 == 1,
            dm = fitsread(dmfitssrc);
        elseif dm1ordm2 == 2,
            dm = fitsread(dmfitssrc, 'image');
        else
            error('dm1 or dm2?');
        end
        return
    end

    %else
    switch class(dm_input)
        case 'double'
            % do nothing
            dm = dm_input;
        
        case 'char'
            dm = fitsread(dm_input);
        
        case 'struct'
            DM_Gain_Map_NMpV = fitsread(PathTranslator(dm_input.gain_fn));
            DM_Gain_Map_NMpV(DM_Gain_Map_NMpV <= 0) = 0;
            dm = (dm_input.dm_displacement/U.NM) ./ DM_Gain_Map_NMpV;
            dm(dm > 100) = 100;
            dm(dm < 0 | isnan(dm)) = 0;
    end
    
end % ValidateDMinput
