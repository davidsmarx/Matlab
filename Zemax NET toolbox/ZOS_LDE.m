classdef ZOS_LDE
    
    properties
        LDE
        Lens = struct;
        %Nsurf % number of surfaces = LDE.NumberOfSurfaces


    end % properties

    methods

        function S = ZOS_LDE(Zapp)
            % S = ZOS_LDE(ZAPP)
            %
            % get ZAPP from Zapp = MATLABZOSConnection(str2double(instance));

            S.LDE = Zapp.PrimarySystem.LDE;
            
            % create the list
            for ns = 0:(S.LDE.NumberOfSurfaces - 1)
                %  LDE row numbers are 0-offset

                %
                ldesurf = S.LDE.GetSurfaceAt(ns);

                % get names of the surfaces and use as variable names
                surflab = char(deblank(ldesurf.Comment.string));
                if isempty(surflab)
                    surflab = ['surf_' num2str(ns)];
                end

                % replace invalid characters with _
                surflab(regexp(surflab,'\W')) = '_';
                % cannot start with a number
                if regexp(surflab(1),'\d'), surflab = ['s_' surflab]; end
                
                % check for identical surface labels
                if isfield(S.Lens,surflab)
                    surflab = [surflab '_' num2str(ns)];
                    %disp(['Warning: identical surface name, changing to: ' surflab]);
                end
   
                % add to the Lens struct
                try
                    S.Lens.(surflab) = ldesurf;
                    %eval(['Lens.' surflab '= surfst;']);
                catch ME
                    disp(ME.stack);
                    disp(ME.message);
                end

            end % for each surface

        end % init ZOS_LDE

        function [Rot, Tran] = GlobalMatrix(S, surfnumber)
            % Get global matrix from LDE, and 
            % parse into Rotation matrix Rot
            % and translation vector Tran in (m)
            %
            % e.g. raypos_global = Rot * raypos(:) + Tran;
        
            MM = 1e-3;

            [success, Rg(1,1), Rg(1,2), Rg(1,3),...
                Rg(2,1), Rg(2,2), Rg(2,3),...
                Rg(3,1), Rg(3,2), Rg(3,3),...
                Rg(4,1), Rg(4,2), Rg(4,3)] = S.LDE.GetGlobalMatrix(surfnumber);
            
            Rot = Rg(1:3,1:3);
            Tran = Rg(4,:).'*MM; % we're working in mks

            if ~success
                error(['GetGlobalMatrix() Failed: ' num2str(success)]);
            end

        end % GlobalMatrix

    end % methods

end % classdef

