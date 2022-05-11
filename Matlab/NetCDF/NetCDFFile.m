classdef NetCDFFile < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        path
        ncid

        dimensions
        variables
        attributes
    end

    methods
        function self = NetCDFFile(path,overwriteExisting)
            self.path = path;
            shouldOverwrite = 0;
            if nargin == 2 && strcmp(overwriteExisting,'OVERWRITE_EXISTING')
                shouldOverwrite = 1;
            end
            if isfile(self.path)
                if shouldOverwrite == 1
                else
                    self.InitializeFromExistingFile();
                end
            end
        end

        function InitializeFromExistingFile(self)
            self.open();

            [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(self.ncid);

            self.dimensions = NetCDFDimension.empty(ndims,0);
            for iDim=0:(ndims-1)
                [dimname,dimlen] = netcdf.inqDim(self.ncid,iDim);
                self.dimensions(iDim+1) = NetCDFDimension(dimname,dimlen,iDim);
                if iDim == unlimdimid
                    self.dimensions(iDim+1).isMutable = 1;
                end
            end
            self.dimensions = self.dimensions.';

            self.variables = NetCDFVariable.empty(nvars,0);
            for iVar=0:(nvars-1)
                [varname, xtype, dimids, numatts] = netcdf.inqVar(self.ncid,iVar);
                variableAttributes = containers.Map;
                for iAtt=0:(numatts-1)
                    gattname = netcdf.inqAttName(self.ncid,iVar,iAtt);
                    variableAttributes(gattname) = netcdf.getAtt(self.ncid,iVar,gattname);
                end
                self.variables(iVar+1) = NetCDFVariable(varname,self.dimensionsForDimIDs(dimids),variableAttributes,self.typeStringForTypeID(xtype),iVar);
            end
            self.variables = self.variables.';

            self.attributes = containers.Map;
            for iAtt=0:(ngatts-1)
                gattname = netcdf.inqAttName(self.ncid,netcdf.getConstant('NC_GLOBAL'),iAtt);
                self.attributes(gattname) = netcdf.getAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'),gattname);
            end
        end

        function [dimension,variable] = addCoordinateDimension(self,data,name,properties)
            
        end

        function variable = addVariable(self,data,name,dims,properties)

        end

        function dims = dimensionsForDimIDs(self,dimids)
            dims = NetCDFDimension.empty(length(dimids),0);
            for iDim=1:length(dimids)
                dims(iDim) = self.dimensions(dimids(iDim)+1);
            end
        end

        function val = typeStringForTypeID(self,type)
            types = {'NC_DOUBLE','NC_FLOAT','NC_INT64','NC_UINT64','NC_INT','NC_UINT','NC_SHORT','NC_USHORT','NC_BYTE','NC_UBYTE','NC_CHAR','NC_CHAR'};
            val = nan;
            for i=1:length(types)
                if netcdf.getConstant(types{i}) == type
                    val = types{i};
                    return
                end
            end
        end

        function self = sync(self)
            netcdf.sync(self.ncid);
        end

        function self = open(self)
            self.ncid = netcdf.open(self.path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
        end

        function self = close(self)
            netcdf.close(self.ncid);
            self.ncid = [];
        end
    end
end