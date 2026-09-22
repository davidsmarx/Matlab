classdef ZOS_HuygensPsf

    properties

        HuygensPSF
        SettingsCfg_fn
        FieldNumber

    end % properties

    properties (Constant)
        MM = 1e-3;
        UM = 1e-6;
    end

    methods

        function S = ZOS_HuygensPsf(Zapp, varargin)

            % the analysis object
            S.HuygensPSF = Zapp.PrimarySystem.Analyses.New_HuygensPsf();

            % options
            S.SettingsCfg_fn = CheckOption('settingscfg_fn', [], varargin{:}); % must be full path

            % 
            if ~isempty(S.SettingsCfg_fn)
                bsuccess = S.HuygensPSF.Settings.LoadFrom(S.SettingsCfg_fn);
                if ~bsuccess,
                    error(['error loading config from ' S.SettingsCfg_fn '\ncheck fn includes full path.']);
                end
            end

            S.FieldNumber = S.HuygensPSF.Settings.Field.GetFieldNumber;

        end % init

        function fieldnum = get.FieldNumber(S)
            %fieldnum = S.HuygensPSF.Settings.Field.GetFieldNumber;
            fieldnum = S.FieldNumber;

        end % get FieldNumber

        function S = set.FieldNumber(S, fieldnum)
            S.HuygensPSF.Settings.Field.SetFieldNumber(fieldnum);
            S.FieldNumber = S.HuygensPSF.Settings.Field.GetFieldNumber;
        end % get FieldNumber

        function [psf, x, y, headerinfo] = Calculate(S)

            S.HuygensPSF.ApplyAndWaitForCompletion();

            % for convenience
            data = S.HuygensPSF.Results.DataGrids(1);

            psf = double(data.ValueData.Data);
            nx = double(data.Nx);
            ny = double(data.Ny);
            x = data.MinX + (0:(nx-1))*data.Dx*S.UM;
            y = data.MinY + (0:(ny-1))*data.Dy*S.UM;

            % textfile to get headerinfo
            fn_text = fullfile(pwd, 'huygenspsf.txt');
            S.HuygensPSF.Results.GetTextFile(fn_text);
            [field, headerinfo] = z_GetPSFHuygens('textfilename', fn_text);
            % psf = flipud(field);

        end % Calculate

    end % methods

end % classdef


%%%%%% example
% huyg = ZOS_HuygensPsf();
% methodsview(huyg)
% 
% % field point
% huyg.Settings.Field
% huyg.Settings.Field.SetFieldNumber(4)
% huyg.Settings.Field.GetFieldNumber
% 
% % calculate
% huyg.ApplyAndWaitForCompletion();
% huyg.Results
% 
% 
%  huyg.Results.DataGrids(1)
% 
% ans = 
% 
%   AR_DataGrid with properties:
% 
%     Description: [1×1 System.String]
%              Dx: 2.6178
%              Dy: 2.6178
%            MinX: -83.7695
%            MinY: -83.7695
%          XLabel: []
%          YLabel: []
%      ValueLabel: []
%       ValueData: [1×1 ZemaxUI.ZOSAPI.SerializedMatrixData]
%          Values: [1×1 System.Double[,]]
%              Nx: 64
%              Ny: 64
%
% % get the data array
% aaa = double(huyg.Results.DataGrids(1).ValueData.Data);
% 
