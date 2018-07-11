classdef CDataHandle < handle
    % this class is for the purpose of creating variables that behave like
    % pointers. For example, multiple objects can set and read the value.
    properties
        dataval;
        
    end % properties
    
    methods
    
        function S = CDataHandle
            S.dataval = [];
        
        end
        
    end % methods
    
end % classdef

