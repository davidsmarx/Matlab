classdef ZOS_Fields

    properties

        Cfields
        units

    end % properties

    properties (Constant)
        MM = 1e-3;
        P = pi/180; % degrees
        PARCMIN = pi/180/60;
        PARCSEC = pi/180/60/60;
    end

    methods

        function S = ZOS_Fields(Zapp)

            S.Cfields = Zapp.PrimarySystem.SystemData.Fields;
            fieldtype = S.Cfields.GetFieldType;

            % units so input/output are mks
            switch char(fieldtype)
                case 'Angle'
                    S.units = S.P;
                case 'ObjectHeight'
                    S.units = S.MM;
                otherwise
                    error(['Field Type: ' char(fieldtype)]);
            end

        end % init

        function S = set_Field(S, ifield, val_xy)
            % val_xy = mks units

            S.Cfields.GetField(ifield).X = val_xy(1)./S.units;
            S.Cfields.GetField(ifield).Y = val_xy(2)./S.units;

        end % set_Field

        function val_xy = get_Field(S, ifield)
            % val_xy = mks units

            val_xy = [
                S.Cfields.GetField(ifield).X*S.units
                S.Cfields.GetField(ifield).Y*S.units
                ];
            
        end % set_Field
        
        function [max_field, i_field] = MaxField(S)
            % [max_field, i_field] = MaxField(S)
            % max_field (radians or m)

            max_field = -1E9;

            for i_field = 1:S.Cfields.NumberOfFields

                xy_field = S.get_Field(i_field);
                radial_field= hypot(xy_field(1), xy_field(2));

                if (radial_field> max_field)

                    max_field = radial_field;
                    max_field_num=i_field;

                end

            end % for each field

            fprintf('max field norm radius: %d, %f\n', max_field_num, max_field);

        end % MaxField

    end % methods

end % classdef 
       