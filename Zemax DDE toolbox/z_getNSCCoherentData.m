function DetData = z_getNSCCoherentData(zchan, nsurf, nobject, pixel, datatype)
% DetData = z_getNSCCoherentData(zchan, nsurf, nobject, pixel, datatype)
%
% nsurf = surface number for NSC in Lens Data Editor
% nobject = object number in NSC
% npixel = pixel # on the detector, or = 0 for sum of all pixels
% datatype = '0' => real, '1' => imaginary, '2' => amplitude, '3' => power

%     surf = ',1';
%     object = sprintf(',%d', sDetectorParms.DetectorObject);
%     pixel = ',0'; % return sum of all pixels for the detector
%     data = ',3'; % return power
% NSCCoherentData,surface,object,pixel,data
    cmdstr = ['NSCCoherentData,' num2str(nsurf)...
        ',' num2str(nobject)...
        ',' num2str(pixel)...
        ',' datatype];
    
    stmp = ddereq(zchan,cmdstr,[1 1]);
    
    DetData = str2double(stmp);

return
