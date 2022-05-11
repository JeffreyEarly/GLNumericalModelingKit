classdef NetCDFVariable
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        varID
        name
        dimensions
        attributes
        isComplex
    end

    methods
        function self = NetCDFVariable(name,dimensions,attributes,varID)
            self.name = name;
            self.dimensions = dimensions;
            self.attributes = attributes;
            self.varID = varID;
        end

    end
end