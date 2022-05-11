classdef NetCDFFile < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        path
        ncid

        dimensions
        variables
        attributes

        dimensionWithName
        variableWithName

        complexVariables
        complexVariableWithName
    end

    properties (Constant)
        GLNetCDFSchemaVersionKey = "GLNetCDFSchemaVersion";

        % attributes for coordinate variables (dimensions)
        GLNetCDFSchemaIsCoordinateVariableKey = "isCoordinateVariable";
        GLNetCDFSchemaIsPeridiocKey = "isPeriodic";
        GLNetCDFSchemaMutableKey = "isMutable";
        GLNetCDFSchemaDomainMinimumKey = "domainMin";
        GLNetCDFSchemaBasisFunctionKey = "basisFunction";
        GLNetCDFSchemaDomainLengthKey = "domainLength";
        GLNetCDFSchemaIsEvenlySampledKey = "isEvenlySampled";
        GLNetCDFSchemaSampleIntervalKey = "sampleInterval";
        GLNetCDFSchemaGridTypeKey = "gridType";
        GLNetCDFSchemaIsFrequencyDomainKey = "isFrequencyDomain";

        % attributes for variables and dimensions
        GLNetCDFSchemaUnitsKey = "units";

        % attributes for variables
        GLNetCDFSchemaIsComplexKey = "isComplex";
        GLNetCDFSchemaProperNameKey = "properName";
        GLNetCDFSchemaIsRealPartKey = "isRealPart";
        GLNetCDFSchemaIsImaginaryPartKey = "isImaginaryPart";
        
        GLNetCDFSchemaUniqueVariableIDKey = "uniqueVariableID";
    end

    methods
        function self = NetCDFFile(path,overwriteExisting)
            self.path = path;
            shouldOverwrite = 0;
            self.dimensions = NetCDFDimension.empty(0,0);
            self.dimensionWithName = containers.Map;
            self.variables = NetCDFVariable.empty(0,0);
            self.variableWithName = containers.Map;
            self.attributes = containers.Map;
            self.complexVariables = NetCDFComplexVariable.empty(0,0);
            self.complexVariableWithName = containers.Map;

            if nargin == 2 && strcmp(overwriteExisting,'OVERWRITE_EXISTING')
                shouldOverwrite = 1;
            end
            if isfile(self.path)
                if shouldOverwrite == 1
                    delete(self.path);
                    self.CreateNewFile();
                else
                    self.InitializeFromExistingFile();
                end
            else
                self.CreateNewFile();
            end
        end

        function CreateNewFile(self)
            self.ncid = netcdf.create(self.path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
            netcdf.endDef(self.ncid);
        end

        function InitializeFromExistingFile(self)
            self.open();

            [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(self.ncid);

            for iDim=0:(ndims-1)
                [dimname,dimlen] = netcdf.inqDim(self.ncid,iDim);
                self.dimensions(iDim+1) = NetCDFDimension(dimname,dimlen,iDim);
                self.dimensionWithName(dimname) = self.dimensions(iDim+1);
                if iDim == unlimdimid
                    self.dimensions(iDim+1).isMutable = 1;
                end
            end
            self.dimensions = self.dimensions.';

            for iVar=0:(nvars-1)
                [varname, xtype, dimids, numatts] = netcdf.inqVar(self.ncid,iVar);
                variableAttributes = containers.Map;
                for iAtt=0:(numatts-1)
                    gattname = netcdf.inqAttName(self.ncid,iVar,iAtt);
                    variableAttributes(gattname) = netcdf.getAtt(self.ncid,iVar,gattname);
                end
                self.variables(iVar+1) = NetCDFVariable(varname,self.dimensionsForDimIDs(dimids),variableAttributes,self.typeStringForTypeID(xtype),iVar);
                self.variableWithName(varname) = self.variables(iVar+1);
            end
            self.variables = self.variables.';
 
            for iAtt=0:(ngatts-1)
                gattname = netcdf.inqAttName(self.ncid,netcdf.getConstant('NC_GLOBAL'),iAtt);
                self.attributes(gattname) = netcdf.getAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'),gattname);
            end
        end

        function [dimension,variable] = addMutableDimension(self,name,properties)
            [dimension,variable] = self.addDimension([],name,properties,1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Dimensions
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [dimension,variable] = addDimension(self,data,name,properties,isMutable)
            
            if nargin < 5
                isMutable = 0;
            end
            if isKey(self.dimensionWithName,name) || isKey(self.variableWithName,name)
                error('A dimension with that name already exists.');
            end
            if (isMutable == 0 && isempty(data)) || isempty(name)
                error('You must specify a name and data');
            end
            if isMutable == 1
                n = netcdf.getConstant('NC_UNLIMITED');
            else
                n = length(data);
            end
            
            ncType = self.netCDFTypeForData(data);
            netcdf.reDef(self.ncid);
            dimID = netcdf.defDim(self.ncid, name, n);
            varID = netcdf.defVar(self.ncid, name, ncType, dimID);
            if ~isempty(properties)
                keyNames = keys(properties);
                for iKey = 1:length(keyNames)
                    netcdf.putAtt(self.ncid,varID, keyNames{iKey}, properties(keyNames{iKey}));
                end
            end
            netcdf.endDef(self.ncid);
            if isMutable == 0 && ~isempty(data)
                netcdf.putVar(self.ncid, varID, data);
            end
            
            dimension = NetCDFDimension(name,n,dimID);
            variable = NetCDFVariable(name,dimension,properties,ncType,varID);

            % Can this fail? I don't think it should.
            self.dimensions(dimID+1) = dimension;
            self.dimensionWithName(name) = dimension;
            self.variables(varID+1) = variable;
            self.variableWithName(name) = variable;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function complexVariable = initComplexVariable(self,name,dimNames,ncType,properties)
            if isKey(self.complexVariableWithName,name)
                error('A variable with that name already exists.');
            end
            if nargin < 5 || isempty(properties)
                properties = containers.Map;
            end
            properties(self.GLNetCDFSchemaIsComplexKey) = 1;

            properties(self.GLNetCDFSchemaIsRealPartKey) = 1;
            properties(self.GLNetCDFSchemaIsImaginaryPartKey) = 0;
            realVar = self.initVariable(strcat(name,"_realp"),dimNames,ncType,properties);

            properties(self.GLNetCDFSchemaIsRealPartKey) = 0;
            properties(self.GLNetCDFSchemaIsImaginaryPartKey) = 1;
            imagVar = self.initVariable(strcat(name,"_imagp"),dimNames,ncType,properties);
            
            complexVariable = NetCDFComplexVariable(name,realVar,imagVar);
            self.complexVariables(end+1) = complexVariable;
            self.complexVariableWithName(name) = complexVariable;
        end

        function variable = initVariable(self,name,dimNames,ncType,properties)
            if isKey(self.variableWithName,name)
                error('A variable with that name already exists.');
            end
            dimIDs = nan(length(dimNames),1);
            for iDim=1:length(dimNames)
                if ~isKey(self.dimensionWithName,dimNames{iDim})
                    error('Unable to find a dimension with the name %s',dimNames{iDim});
                end
                dimIDs(iDim) = self.dimensionWithName(dimNames{iDim}).dimID;
            end

            netcdf.reDef(self.ncid);
            varID = netcdf.defVar(self.ncid, name, ncType, dimIDs);
            if ~isempty(properties)
                keyNames = keys(properties);
                for iKey = 1:length(keyNames)
                    netcdf.putAtt(self.ncid,varID, keyNames{iKey}, properties(keyNames{iKey}));
                end
            end
            netcdf.endDef(self.ncid);

            variable = NetCDFVariable(name,self.dimensionsForDimIDs(dimIDs),properties,ncType,varID);
            self.variables(varID+1) = variable;
            self.variableWithName(name) = variable;
        end

        function setVariable(data,name)
            if isKey(self.complexVariableWithName,name)
                complexVariable = self.complexVariableWithName(name);
                variable = complexVariable.realVar;
                for iDim=1:length(variable.dimensions)
                    if size(data,iDim) ~= variable.dimensions(iDim).nPoints
                        error('Incorrect dimension size: dimension %d of the data is length %d, but the dimension %s has length %d.',iDim,size(data,iDim),variable.dimensions(iDim).name,variable.dimensions(iDim).nPoints);
                    end
                end
                netcdf.putVar(self.ncid, complexVariable.realVar.varID, real(data));
                netcdf.putVar(self.ncid, complexVariable.imagVar.varID, imag(data));
            else
                variable = self.variableWithName(name);
                for iDim=1:length(variable.dimensions)
                    if size(data,iDim) ~= variable.dimensions(iDim).nPoints
                        error('Incorrect dimension size: dimension %d of the data is length %d, but the dimension %s has length %d.',iDim,size(data,iDim),variable.dimensions(iDim).name,variable.dimensions(iDim).nPoints);
                    end
                end
                netcdf.putVar(self.ncid, variable.varID, data);
            end

        end

        function variable = addVariable(self,data,name,dimNames,properties)
            if isreal(data)
                ncType = self.netCDFTypeForData(data);
                variable = self.initVariable(name,dimNames,ncType,properties);
                self.setVariable(data,name);
            else
                ncType = self.netCDFTypeForData(data);
                variable = self.complexVariable(name,dimNames,ncType,properties);
                self.setVariable(data,name);
            end
        end

        function concatenateVariableAlongDimension(self,data,name,dimName,index)
            if isKey(self.complexVariableWithName,name)
                complexVariable = self.complexVariableWithName(name);
                variable = complexVariable.realVar;
                start = zeros(1,length(variable.dimensions));
                count = zeros(1,length(variable.dimensions));
                for iDim=1:length(variable.dimensions)
                    if strcmp(variable.dimensions(iDim).name,dimName)
                        start(iDim) = index-1;
                        count(iDim) = 1;
                    else
                        start(iDim) = 0;
                        count(iDim) = variable.dimensions(iDim).nPoints;
                    end
                end
                netcdf.putVar(self.ncid, complexVariable.realVar.varID, start, count, real(data));
                netcdf.putVar(self.ncid, complexVariable.imagVar.varID, start, count, imag(data));
            elseif isKey(self.variableWithName,name)
                variable = self.variableWithName(name);
                start = zeros(1,length(variable.dimensions));
                count = zeros(1,length(variable.dimensions));
                for iDim=1:length(variable.dimensions)
                    if strcmp(variable.dimensions(iDim).name,dimName)
                        start(iDim) = index-1;
                        count(iDim) = 1;
                    else
                        start(iDim) = 0;
                        count(iDim) = variable.dimensions(iDim).nPoints;
                    end
                end
                netcdf.putVar(self.ncid, variable.varID, start, count, data);
            else
                error('Unable to find a variable with the name %s',name);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Utilities
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function dims = dimensionsForDimIDs(self,dimids)
            dims = NetCDFDimension.empty(length(dimids),0);
            for iDim=1:length(dimids)
                dims(iDim) = self.dimensions(dimids(iDim)+1);
            end
        end

        function val = netCDFTypeForData(self,data)
            keys = {'double','single','int64','uint64','int32','uint32','int16','uint16','int8','uint8','char','string'};
            values = {'NC_DOUBLE','NC_FLOAT','NC_INT64','NC_UINT64','NC_INT','NC_UINT','NC_SHORT','NC_USHORT','NC_BYTE','NC_UBYTE','NC_CHAR','NC_CHAR'};
            map = containers.Map(keys, values);
            if ~isKey(map,class(data))
                error('unknown data type');
            end
            val = map(class(data));
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