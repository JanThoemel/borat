clc; clear all; close all;

set(groot,'DefaultFigureColormap',gray)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NAME='BlackOut RAy-Tracer';NAMESHORT='BORAT ';VERSION=3.0;

% ExoMars UHF frequency (TBC)

f=400.0e6;                  

x_0=1.252+0.001;
y_0=0;
z_0=1.066+0.001;

cut_off_limit=0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha1=0;
alpha2=360;

n_alpha=2; %n_point per ring, longitudinal directions

n_theta=2; %n_ring from azimuth

semi_aperture=70; %semi-aperture of antenna


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%fix box for reflection

max_reflections=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Input files

%tecplot dat file with format x-y-z-rho_electrons-T-Ri-dRi/dx-dRi/dy-dRi/dz

% input_method = 'binary tecplot file';
% flowfieldfilename='CFDsolutions/ExoMars_3D_62s.szplt';
% wallfilename='3Dmeshes/ExoMars_capsule.szplt';

input_method = 'ascii dat file';
wallfilename='3Dmeshes/ExoMars_capsule.dat';
flowfieldfilename='CFDsolutions/ExoMars_3D_62s.dat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Output file

solution_file='ExoMars solution';

figure_file = 'RayTracing solution';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Read solution file

switch input_method
    
    % ascii
    case 'ascii dat file'
        z=0;
        domain.nova=9; % x-y-z-rho_electrons-T-Ri-dRi/dx-dRi/dy-dRi/dz
        id=fopen(flowfieldfilename);
        %read titleline, to be discarded
        dimline=fgetl(id);
        %read variables' name line, to be discarded, but number of variables to be
        %determined by number of quotation marks
        
        while ~feof(id)
            z=z+1;
            
            while 1
                dimline=fgetl(id);
                if contains(dimline,'ZONE')
                    break
                end
            end
            
            while 1
                dimline=fgetl(id);
                if contains(dimline,'Nodes')
                    
                    Npos=strfind(dimline,'Nodes=');
                    Epos=strfind(dimline,'Elements=');
                    Fpos=strfind(dimline,'ZONETYPE=');
                    
                    domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+6:Epos-3));
                    domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+9:Fpos-3));
                    
                    break
                end
            end
            
            while 1
                dimline=fgetl(id);
                if contains(dimline,'DT=')
                    break
                end
            end
            
            domain.( strcat('zone',num2str(z)) ).variables=zeros(  domain.nova  ,  domain.( strcat('zone',num2str(z)) ).N);
            
            
            
            for i=1:domain.( strcat('zone',num2str(z)) ).N
                line=fgetl(id);
                domain.( strcat('zone',num2str(z)) ).variables(1:domain.nova,i)=sscanf(line,'  %g');
            end
            
            
            domain.( strcat('zone',num2str(z)) ).NodeMap=zeros(domain.( strcat('zone',num2str(z)) ).E,4);
            
            for i=1:domain.( strcat('zone',num2str(z)) ).E
                
                line=fgetl(id);
                domain.( strcat('zone',num2str(z)) ).NodeMap(i,:)=sscanf(line,'  %g');
                
            end
            
            
        end % while eof
      
        fclose(id);
        domain.nozones=z;
        
        
        domain.( strcat('zone',num2str(z)) ).DT=triangulation(domain.( strcat('zone',num2str(z)) ).NodeMap,domain.( strcat('zone',num2str(z)) ).variables(1:3,:)');
        
        [domain.( strcat('zone',num2str(z)) ).BoundaryFaces,domain.( strcat('zone',num2str(z)) ).BoundaryPoints ]=...
            freeBoundary(domain.( strcat('zone',num2str(z)) ).DT);
        
        % wall file
        z=0;
        domain.WallNova=3; %just 3D points
        id=fopen(wallfilename);
        %read titleline, to be discarded
        dimline=fgetl(id);
        %read variables' name line, to be discarded, but number of variables to be
        %determined by number of quotation marks
        
        while ~feof(id)
            z=z+1;
            
            while 1
                dimline=fgetl(id);
                if contains(dimline,'ZONE')
                    break
                end
            end
            
            while 1
                dimline=fgetl(id);
                if contains(dimline,'Nodes')
                    
                    Npos=strfind(dimline,'Nodes=');
                    Epos=strfind(dimline,'Elements=');
                    Fpos=strfind(dimline,'ZONETYPE=');
                    
                    domain.( strcat('zone',num2str(z)) ).Wall.N=str2double(dimline(Npos+6:Epos-3));
                    domain.( strcat('zone',num2str(z)) ).Wall.E=str2double(dimline(Epos+9:Fpos-3));
                    
                    break
                end
            end
            
            while 1
                dimline=fgetl(id);
                if contains(dimline,'DT=')
                    break
                end
            end
            
            domain.( strcat('zone',num2str(z)) ).Wall.variables=zeros(  domain.WallNova  ,  domain.( strcat('zone',num2str(z)) ).Wall.N);
            
            
            
            for i=1:domain.( strcat('zone',num2str(z)) ).Wall.N
                line=fgetl(id);
                domain.( strcat('zone',num2str(z)) ).Wall.variables(1:domain.WallNova,i)=sscanf(line,'  %g');
            end
            
            
            domain.( strcat('zone',num2str(z)) ).Wall.NodeMap=zeros(domain.( strcat('zone',num2str(z)) ).Wall.E,3);
            
            for i=1:domain.( strcat('zone',num2str(z)) ).Wall.E
                
                line=fgetl(id);
                domain.( strcat('zone',num2str(z)) ).Wall.NodeMap(i,:)=sscanf(line,'  %g');
                
            end
            
            
            
        end % while eof
        fclose(id);
        domain.nozones=z;
        
        
        domain.(strcat('zone',num2str(z))).Wall.BoundaryDT=triangulation(domain.(strcat('zone',num2str(z))).Wall.NodeMap,domain.(strcat('zone',num2str(z))).Wall.variables(1:3,:)');
        domain.( strcat('zone',num2str(z)) ).Wall.CenterPoint=incenter(domain.( strcat('zone',num2str(z)) ).Wall.BoundaryDT);
        
        
    case 'binary tecplot file'
        
        tecplot_home = '/Applications/Tecplot 360 EX 2020 R2/Tecplot 360 EX 2020 R2.app/Contents';
        tecio_path = strcat(tecplot_home, '/Frameworks/libtecio.dylib');
        tecio_header_path = '/TecIO/TECIO.h'; %Uses simplified TECIO.h included in package.
        
        if ~libisloaded('tecio')
            [notfound, warnings] = loadlibrary(tecio_path, tecio_header_path,...
                'alias', 'tecio');
        end
        
        % Open solution szplt file
        
        % Pass NULL by [] to get the filehandle
        [isok, ~, handle] = calllib('tecio', 'tecFileReaderOpen', flowfieldfilename, []);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        dataSetTitle = libpointer('stringPtrPtr', cell(1,1));
        [isok, ~, dataSetTitle] = calllib('tecio', 'tecDataSetGetTitle', handle,  dataSetTitle);
        
        domain.dataTitle=char(dataSetTitle);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fileType=0;
        [isok, ~, fileType] = calllib('tecio', 'tecFileGetType', handle, fileType);
        
        domain.fileType=fileType;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        numZones = 0; % Defining variable
        [isok, ~, numZones] = calllib('tecio', 'tecDataSetGetNumZones', handle,  numZones);
        domain.nozones=numZones;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        numVars = 0;
        [isok, ~, numVars] = calllib('tecio', 'tecDataSetGetNumVars', handle,  numVars);
        domain.nova=numVars;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        domain.VarNamesList=[];
        
        for var=1:numVars
            
            VarNames = libpointer('stringPtrPtr', cell(1,1));
            
            [isok, ~, VarNames] = calllib('tecio', 'tecVarGetName', handle, var,  VarNames);
            
            domain.VarNamesList = strcat(domain.VarNamesList,char(VarNames));
            if var<numVars
                domain.VarNamesList = strcat(domain.VarNamesList,',');
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        domain.unstructured=1; %to be checked
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for z = 1:numZones
            %%Get information on what the Zone Type is.
            
            zoneType= 0;
            [isok, ~, zoneType] = calllib('tecio', 'tecZoneGetType', handle, z,zoneType);
            domain.( strcat('zone',num2str(z)) ).zoneType=zoneType;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            zone_title = libpointer('stringPtrPtr', cell(1,1));
            [isok, ~, zone_title] = calllib('tecio', 'tecZoneGetTitle', handle, z,  zone_title);
            
            domain.( strcat('zone',num2str(z)) ).zoneTitle = char(zone_title);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            iMax = 0; jMax = 0; kMax = 0;
            [isok, ~, iMax, jMax, kMax] = calllib('tecio', 'tecZoneGetIJK', handle, z, iMax, jMax, kMax);
            
            domain.( strcat('zone',num2str(z)) ).iMax=iMax;
            domain.( strcat('zone',num2str(z)) ).jMax=jMax;
            domain.( strcat('zone',num2str(z)) ).kMax=kMax;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            varTypes = zeros(numVars,1);
            passiveVarList = zeros(numVars,1);
            valueLocation = zeros(numVars,1);
            shareVarFromZone = zeros(numVars,1);
            
            for var = 1:numVars
                
                [isok, ~, varTypes(var,1)] = calllib('tecio', 'tecZoneVarGetType',...
                    handle, z,var, varTypes(var,1));
                
                [isok, ~, passiveVarList(var,1)] = calllib('tecio', 'tecZoneVarGetSharedZone',...
                    handle, z,var, passiveVarList(var,1));
                
                [isok, ~, valueLocation(var,1)] = calllib('tecio', 'tecZoneVarGetValueLocation',...
                    handle, z,var, valueLocation(var,1));
                
                [isok, ~, shareVarFromZone(var,1)] = calllib('tecio', 'tecZoneVarIsPassive',...
                    handle, z,var, shareVarFromZone(var,1));
                
            end
            
            domain.( strcat('zone',num2str(z)) ).varTypes=varTypes;
            domain.( strcat('zone',num2str(z)) ).passiveVarList=passiveVarList;
            domain.( strcat('zone',num2str(z)) ).valueLocation=valueLocation;
            domain.( strcat('zone',num2str(z)) ).shareVarFromZone=shareVarFromZone;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            shareConnectivityFromZone =0;
            faceNeighborMode =0;
            numFaceConnections =0;
            
            
            [isok, ~, shareConnectivityFromZone] = calllib('tecio', 'tecZoneConnectivityGetSharedZone',...
                handle, z, shareConnectivityFromZone);
            
            [isok, ~, faceNeighborMode] = calllib('tecio', 'tecZoneFaceNbrGetMode',...
                handle, z, faceNeighborMode);
            
            [isok, ~, numFaceConnections] = calllib('tecio', 'tecZoneFaceNbrGetNumConnections',...
                handle, z, numFaceConnections);
            
            
            domain.( strcat('zone',num2str(z)) ).shareConnectivityFromZone=shareConnectivityFromZone;
            domain.( strcat('zone',num2str(z)) ).faceNeighborMode=faceNeighborMode;
            domain.( strcat('zone',num2str(z)) ).numFaceConnections=numFaceConnections;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            for var = 1:numVars
                
                numValues=0;
                [isok, ~, numValues] = calllib('tecio', 'tecZoneVarGetNumValues',...
                    handle, z,var, numValues);
                
            end
            
            domain.( strcat('zone',num2str(z)) ).N=numValues;
            domain.( strcat('zone',num2str(z)) ).variables=zeros(  domain.nova  ,  domain.( strcat('zone',num2str(z)) ).N);
            
            for var = 1:numVars
                
                if varTypes(var,1)==1
                    [isok, ~, domain.( strcat('zone',num2str(z)) ).variables(var,:)] = calllib('tecio', 'tecZoneVarGetFloatValues', ...
                        handle, z,var, 1,numValues,domain.( strcat('zone',num2str(z)) ).variables(var,:));
                    
                else
                    [isok, ~, domain.( strcat('zone',num2str(z)) ).variables(var,:)] = calllib('tecio', 'tecZoneVarGetDoubleValues', ...
                        handle, z,var, 1,numValues,domain.( strcat('zone',num2str(z)) ).variables(var,:));
                end
                
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            if (zoneType ~= 0 && shareConnectivityFromZone == 0)
                
                
                numValues=0;
                [isok, ~, numValues] = calllib('tecio', 'tecZoneNodeMapGetNumValues', ...
                    handle, z,jMax,numValues);
                
                is64Bit= 0;
                [isok, ~, is64Bit] = calllib('tecio', 'tecZoneNodeMapIs64Bit', handle, z, is64Bit);
                
                domain.( strcat('zone',num2str(z)) ).E=numValues;
                domain.( strcat('zone',num2str(z)) ).is64Bit=is64Bit;
                
                if is64Bit
                    
                    nodeMap=zeros(numValues,1);
                    
                    [isok, ~, nodeMap] = calllib('tecio', 'tecZoneNodeMapGet64',...
                        handle, z,1, jMax, nodeMap(:,1));
                    
                    
                    domain.( strcat('zone',num2str(z)) ).originalNodeMap=nodeMap;
                    
                else
                    
                    nodeMap=zeros(numValues,1);
                    [isok, ~, nodeMap] = calllib('tecio', 'tecZoneNodeMapGet',...
                        handle, z,1, jMax, nodeMap(:,1));
                    
                    domain.( strcat('zone',num2str(z)) ).originalNodeMap=nodeMap;
                    
                    
                    if domain.( strcat('zone',num2str(z)) ).zoneType==4 % for TETRA mesh
                        
                        domain.( strcat('zone',num2str(z)) ).NodeMap=zeros(jMax,4);
                        
                        domain.( strcat('zone',num2str(z)) ).NodeMap(:,1)=nodeMap(1:4:end);
                        domain.( strcat('zone',num2str(z)) ).NodeMap(:,2)=nodeMap(2:4:end);
                        domain.( strcat('zone',num2str(z)) ).NodeMap(:,3)=nodeMap(3:4:end);
                        domain.( strcat('zone',num2str(z)) ).NodeMap(:,4)=nodeMap(4:4:end);
                        
                    end
                    
                    
                end
                
                
            end
            
            
            
            domain.( strcat('zone',num2str(z)) ).DT=triangulation(domain.( strcat('zone',num2str(z)) ).NodeMap,domain.( strcat('zone',num2str(z)) ).variables(1:3,:)');
            
            [domain.( strcat('zone',num2str(z)) ).BoundaryFaces,domain.( strcat('zone',num2str(z)) ).BoundaryPoints ]=...
                freeBoundary(domain.( strcat('zone',num2str(z)) ).DT);
            
            
            
            
        end %cicle in zone
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Import wall spacecraft
        
        
        % Pass NULL by [] to get the filehandle
        [isok, ~, handle] = calllib('tecio', 'tecFileReaderOpen', wallfilename, []);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        dataSetTitle = libpointer('stringPtrPtr', cell(1,1));
        [isok, ~, dataSetTitle] = calllib('tecio', 'tecDataSetGetTitle', handle,  dataSetTitle);
        
        domain.dataTitle=char(dataSetTitle);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fileType=0;
        [isok, ~, fileType] = calllib('tecio', 'tecFileGetType', handle, fileType);
        
        domain.fileType=fileType;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        numZones = 0; % Defining variable
        [isok, ~, numZones] = calllib('tecio', 'tecDataSetGetNumZones', handle,  numZones);
        domain.nozones=numZones;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        numVars = 0;
        [isok, ~, numVars] = calllib('tecio', 'tecDataSetGetNumVars', handle,  numVars);
        domain.( strcat('zone',num2str(z)) ).Wall.nova=numVars;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        domain.VarNamesList=[];
        
        for var=1:numVars
            
            VarNames = libpointer('stringPtrPtr', cell(1,1));
            
            [isok, ~, VarNames] = calllib('tecio', 'tecVarGetName', handle, var,  VarNames);
            
            domain.VarNamesList = strcat(domain.VarNamesList,char(VarNames));
            if var<numVars
                domain.VarNamesList = strcat(domain.VarNamesList,',');
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        domain.unstructured=1; %to be checked
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for z = 1:numZones
            % Get information on what the Zone Type is.
            
            zoneType= 0;
            [isok, ~, zoneType] = calllib('tecio', 'tecZoneGetType', handle, z,zoneType);
            domain.( strcat('zone',num2str(z)) ).Wall.zoneType=zoneType;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            zone_title = libpointer('stringPtrPtr', cell(1,1));
            [isok, ~, zone_title] = calllib('tecio', 'tecZoneGetTitle', handle, z,  zone_title);
            
            domain.( strcat('zone',num2str(z)) ).Wall.zoneTitle = char(zone_title);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            iMax = 0; jMax = 0; kMax = 0;
            [isok, ~, iMax, jMax, kMax] = calllib('tecio', 'tecZoneGetIJK', handle, z, iMax, jMax, kMax);
            
            domain.( strcat('zone',num2str(z)) ).Wall.iMax=iMax;
            domain.( strcat('zone',num2str(z)) ).Wall.jMax=jMax;
            domain.( strcat('zone',num2str(z)) ).Wall.kMax=kMax;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            varTypes = zeros(numVars,1);
            passiveVarList = zeros(numVars,1);
            valueLocation = zeros(numVars,1);
            shareVarFromZone = zeros(numVars,1);
            
            for var = 1:numVars
                
                [isok, ~, varTypes(var,1)] = calllib('tecio', 'tecZoneVarGetType',...
                    handle, z,var, varTypes(var,1));
                
                [isok, ~, passiveVarList(var,1)] = calllib('tecio', 'tecZoneVarGetSharedZone',...
                    handle, z,var, passiveVarList(var,1));
                
                [isok, ~, valueLocation(var,1)] = calllib('tecio', 'tecZoneVarGetValueLocation',...
                    handle, z,var, valueLocation(var,1));
                
                [isok, ~, shareVarFromZone(var,1)] = calllib('tecio', 'tecZoneVarIsPassive',...
                    handle, z,var, shareVarFromZone(var,1));
                
            end
            
            domain.( strcat('zone',num2str(z)) ).Wall.varTypes=varTypes;
            domain.( strcat('zone',num2str(z)) ).Wall.passiveVarList=passiveVarList;
            domain.( strcat('zone',num2str(z)) ).Wall.valueLocation=valueLocation;
            domain.( strcat('zone',num2str(z)) ).Wall.shareVarFromZone=shareVarFromZone;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            shareConnectivityFromZone =0;
            faceNeighborMode =0;
            numFaceConnections =0;
            
            
            [isok, ~, shareConnectivityFromZone] = calllib('tecio', 'tecZoneConnectivityGetSharedZone',...
                handle, z, shareConnectivityFromZone);
            
            [isok, ~, faceNeighborMode] = calllib('tecio', 'tecZoneFaceNbrGetMode',...
                handle, z, faceNeighborMode);
            
            [isok, ~, numFaceConnections] = calllib('tecio', 'tecZoneFaceNbrGetNumConnections',...
                handle, z, numFaceConnections);
            
            
            domain.( strcat('zone',num2str(z)) ).Wall.shareConnectivityFromZone=shareConnectivityFromZone;
            domain.( strcat('zone',num2str(z)) ).Wall.faceNeighborMode=faceNeighborMode;
            domain.( strcat('zone',num2str(z)) ).Wall.numFaceConnections=numFaceConnections;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            for var = 1:numVars
                
                numValues=0;
                [isok, ~, numValues] = calllib('tecio', 'tecZoneVarGetNumValues',...
                    handle, z,var, numValues);
                
            end
            
            domain.( strcat('zone',num2str(z)) ).Wall.N=numValues;
            domain.( strcat('zone',num2str(z)) ).Wall.variables=zeros(  domain.( strcat('zone',num2str(z)) ).Wall.nova  ,  domain.( strcat('zone',num2str(z)) ).Wall.N);
            
            for var = 1:numVars
                
                if varTypes(var,1)==1
                    [isok, ~, domain.( strcat('zone',num2str(z)) ).Wall.variables(var,:)] = calllib('tecio', 'tecZoneVarGetFloatValues', ...
                        handle, z,var, 1,numValues,domain.( strcat('zone',num2str(z)) ).Wall.variables(var,:));
                    
                else
                    [isok, ~, domain.( strcat('zone',num2str(z)) ).Wall.variables(var,:)] = calllib('tecio', 'tecZoneVarGetDoubleValues', ...
                        handle, z,var, 1,numValues,domain.( strcat('zone',num2str(z)) ).Wall.variables(var,:));
                end
                
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            if (zoneType ~= 0 && shareConnectivityFromZone == 0)
                
                
                numValues=0;
                [isok, ~, numValues] = calllib('tecio', 'tecZoneNodeMapGetNumValues', ...
                    handle, z,jMax,numValues);
                
                is64Bit= 0;
                [isok, ~, is64Bit] = calllib('tecio', 'tecZoneNodeMapIs64Bit', handle, z, is64Bit);
                
                domain.( strcat('zone',num2str(z)) ).Wall.E=numValues;
                domain.( strcat('zone',num2str(z)) ).Wall.is64Bit=is64Bit;
                
                if is64Bit
                    
                    nodeMap=zeros(numValues,1);
                    
                    [isok, ~, nodeMap] = calllib('tecio', 'tecZoneNodeMapGet64',...
                        handle, z,1, jMax, nodeMap(:,1));
                    
                    
                    domain.( strcat('zone',num2str(z)) ).Wall.originalNodeMap=nodeMap;
                    
                else
                    
                    nodeMap=zeros(numValues,1);
                    [isok, ~, nodeMap] = calllib('tecio', 'tecZoneNodeMapGet',...
                        handle, z,1, jMax, nodeMap(:,1));
                    
                    domain.( strcat('zone',num2str(z)) ).Wall.originalNodeMap=nodeMap;
                    
                    
                    domain.( strcat('zone',num2str(z)) ).Wall.NodeMap=zeros(jMax,3);
                    
                    domain.( strcat('zone',num2str(z)) ).Wall.NodeMap(:,1)=nodeMap(1:3:end);
                    domain.( strcat('zone',num2str(z)) ).Wall.NodeMap(:,2)=nodeMap(2:3:end);
                    domain.( strcat('zone',num2str(z)) ).Wall.NodeMap(:,3)=nodeMap(3:3:end);
                    
                end
                
                
            end
            
            
            domain.(strcat('zone',num2str(z))).Wall.BoundaryDT=triangulation(domain.(strcat('zone',num2str(z))).Wall.NodeMap,domain.(strcat('zone',num2str(z))).Wall.variables(1:3,:)');
            domain.( strcat('zone',num2str(z)) ).Wall.CenterPoint=incenter(domain.( strcat('zone',num2str(z)) ).Wall.BoundaryDT);
            
            
            domain.limitx1=min(domain.( strcat('zone',num2str(z)) ).Wall.variables(1,:));
            domain.limitx2=max(domain.( strcat('zone',num2str(z)) ).Wall.variables(1,:));

            domain.limity1=min(domain.( strcat('zone',num2str(z)) ).Wall.variables(2,:));
            domain.limity2=max(domain.( strcat('zone',num2str(z)) ).Wall.variables(2,:));

            domain.limitz1=min(domain.( strcat('zone',num2str(z)) ).Wall.variables(3,:));
            domain.limitz2=max(domain.( strcat('zone',num2str(z)) ).Wall.variables(3,:));
            
        end
        
        
        calllib('tecio', 'tecFileReaderClose', handle) ;
        
        if libisloaded('tecio')
            unloadlibrary('tecio')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    otherwise
        warning('Unexpected import format file.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot Domain

figure;
patch('Vertices',domain.( strcat('zone',num2str(z)) ).Wall.variables(1:3,:)',...
    'Faces',domain.( strcat('zone',num2str(z)) ).Wall.NodeMap,...
    'FaceColor','w','FaceAlpha',1,'EdgeColor','b');
axis tight
view(-51,24)

figure;
patch('Vertices',domain.( strcat('zone',num2str(z)) ).BoundaryPoints,...
    'Faces',domain.( strcat('zone',num2str(z)) ).BoundaryFaces,...
    'FaceColor','w','FaceAlpha',1,'EdgeColor','b');
axis equal
view(-51,24)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Eikonal - Initial conditions


d_alpha=(alpha2-alpha1)/(n_alpha);
d_theta=(semi_aperture)/(n_theta);


for i=1:n_alpha
    
    alpha(i)=  alpha1+d_alpha*(i-1);
    
end

for i=1:n_theta
    
    theta(i)=  90-d_theta*(i);
    
end


total_size=n_theta*n_alpha;

%add slot for absorption
rays_solution_3D=zeros(total_size,1000,8);
rays_reflections=zeros(total_size,1);
rays_refractive_index=zeros(total_size,1000,1);

figure; hold on;
patch('Vertices',domain.( strcat('zone',num2str(z)) ).Wall.variables(1:3,:)',...
    'Faces',domain.( strcat('zone',num2str(z)) ).Wall.NodeMap,...
    'FaceColor','b','FaceAlpha',.7,'EdgeColor','k');
axis tight
axis equal
view(0,0)

ray_counter=1;
for k=1:n_theta
    for i=1:n_alpha
        
        alpha_0=alpha(i);
        %         theta_0=47-(theta(k))*cos(alpha(i)*pi/180)
        theta_0=theta(k);
        
        [n_0,~,~,~] = refractive_interpolation(x_0,y_0,z_0,domain);
       
        chi_x_0=n_0*((cos(theta_0*pi/180))*(cos(alpha_0*pi/180))*cos(47*pi/180)+sin(theta_0*pi/180)*sin(47*pi/180));
        chi_y_0=n_0*(cos(theta_0*pi/180))*(sin(alpha_0*pi/180));
        chi_z_0=n_0*(sin(theta_0*pi/180)*cos(47*pi/180)-(cos(theta_0*pi/180))*(cos(alpha_0*pi/180))*sin(47*pi/180));
        
       
        y0=[x_0,y_0,z_0,chi_x_0,chi_y_0,chi_z_0, n_0,0];
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tstart = 0;
        tfinal =2000;
        
        refine = 4;
        
        
        options = odeset('Events',@(t,y) cutoff(t,y,domain));
                
        options = odeset(options,'Refine',refine); %this doesn'' make it work
        
        options = odeset(options,'InitialStep',1e-2);
        
        options = odeset(options,'MaxStep',1);
        
        tout = tstart;
        yout = y0;
        teout = [];
        yeout = [];
        ieout = [];
        
        
        y_final=[y0(1) y0(2) y0(3)];
        
        flag_out=checkifinsidedomain3D(y_final, domain);
              
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        n_reflections=0;
        
        
        while flag_out
            
            [t,y,te,ye,ie] = ode15s(@(t,y) eikonal(t,y,domain),[tstart tfinal],y0,options);
            %[t,y,te,ye,ie] = ode45(@(t,y) eikonal(t,y,domain),[tstart tfinal],y0,options);
            
            
            % Accumulate output.
            nt = length(t);
            tout = [tout; t(2:nt)];
            yout = [yout; y(2:nt,:)];
            teout = [teout; te];    % Events at tstart are never reported.
            yeout = [yeout; ye];
            ieout = [ieout; ie];
            
            
            
            flag_boundary= checkBoundaryReflection3D([yout(end,1) yout(end,2) yout(end,3)],domain);
            
            %cutoff check
            flag_cutoff=0;
            [n_end, nx_end, ny_end, nz_end] = refractive_interpolation(y(end,1), y(end,2), y(end,3) ,domain);
            
            
            if flag_boundary==1 && n_end>cut_off_limit
                
                n_reflections=n_reflections+1;
                
                if n_reflections>max_reflections
                    flag_out=0;
                    break
                end
                
                P1=[y(end-1,1) y(end-1,2) y(end-1,3)];
                P2=[y(end,1) y(end,2) y(end,3)]; %this is the new starting point
                u=P2-P1;
                
                
                IDboundary=nearestNeighbor(domain.(strcat('zone',num2str(z))).Wall.BoundaryDT,P2);
                N2=vertexNormal(domain.(strcat('zone',num2str(z))).Wall.BoundaryDT,IDboundary);
                
                IDpointCenter=dsearchn(domain.( strcat('zone',num2str(z)) ).Wall.CenterPoint, P2);
                N=faceNormal(domain.(strcat('zone',num2str(z))).Wall.BoundaryDT,IDpointCenter);
                
                if N(1)~=N2(1) && N(2)~=N2(2) && N(3)~=N2(3)
                    N;
                    N2;
                end
                % problems at sharp corners,
                % find other face closer. Is it relevant for finer meshes?
                
                Rr = u - 2*N*(dot(u,N));
                Rr=Rr/sqrt(Rr(1)^2+Rr(2)^2+Rr(3)^2);
                
                
                theta0=asin(Rr(3))*180/pi;
                alpha0=atan(Rr(2)/Rr(1))*180/pi;
                
                if Rr(1)<0
                    alpha0=180+alpha0;
                end
                
                if alpha0<0
                    alpha0=360+alpha0;
                end
                
                
                [n_0,nx_0,ny_0,nz_0] = refractive_interpolation(y(end,1),y(end,2), y(end,3),domain);
                
                y0(4)=n_0*cos(theta0*pi/180)*cos(alpha0*pi/180);
                y0(5)=n_0*cos(theta0*pi/180)*sin(alpha0*pi/180);
                y0(6)=n_0*sin(theta0*pi/180);
                
                y0(1) = y(end,1)+0.01*y0(4);
                y0(2) = y(end,2)+0.01*y0(5);
                y0(3) = y(end,3)+0.01*y0(6);
                
                %  refractive index and path length
                y0(7) = y(end,7);
                y0(8) = y(end,8);
                
                options = odeset(options,'InitialStep',t(nt)-t(nt-refine));
                
                
                tstart = t(nt);
                
                plot3(yout(:,1),yout(:,2), yout(:,3),'-r','Linewidth',1);
                hold on;
                
            elseif n_end<cut_off_limit && checkifinsidedomain3D([yout(end,1) yout(end,2) yout(end,3)], domain)
                
                fprintf('\n \t cut off ');
                flag_out=0;
                
            else
                
                flag_out=0;
                
            end
            
            plot3(yout(:,1),yout(:,2), yout(:,3),'-r','Linewidth',1);
            hold on;
            
        end
        
        fprintf('\n   Ray done N = %g',ray_counter)
        fprintf('   reflections = %g', n_reflections)
        
        %save ray solution
        for istep=1:length(yout)
            rays_solution_3D(ray_counter,istep,:)=yout(istep,:);
            rays_refractive_index(ray_counter,istep,1)=sqrt(yout(istep,4)^2+yout(istep,5)^2+yout(istep,6)^2);
        end
        
        rays_reflections(ray_counter,1)=n_reflections;
        ray_counter=ray_counter+1;
        
    end
    
end


hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savefig(figure_file)
save(solutionfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions eikonal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n,nx,ny,nz] = refractive_interpolation(x,y,z,domain)
    
    isInside=checkifinsidedomain3D([x y z],domain);
    
    if isInside
        
        n=interpolation3D(domain,[x y z],domain.nova-3);
        nx=interpolation3D(domain,[x y z],domain.nova-2);
        ny=interpolation3D(domain,[x y z],domain.nova-1);
        nz=interpolation3D(domain,[x y z],domain.nova);
        
        if n<0.001
            n=0.00;
        end
    else
        
        n=1;
        nx=0;
        ny=0;
        nz=0;
        
    end
        
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function solution = eikonal(t,y,domain)
    
    
    xi=y(1);
    yi=y(2);
    zi=y(3);
    ui=y(4);
    vi=y(5);
    wi=y(6);
    
    if isinf(xi)
        xi;
    end
    
    if isnan(xi)
        xi;
    end
    
    [n, nx, ny,nz] = refractive_interpolation(xi, yi, zi,domain);
    
    if n==0
        n;
        n=0.009;
    end
    
    if n<0.01
        n;
    end
    
    if isinf(n)
        n;
    end
    
    if isnan(n)
        n;
    end
    
    yp = zeros(8,1);
    
    yp(1) = ui / (n^2) ;
    yp(2) = vi / (n^2) ;
    yp(3) = wi / (n^2) ;
    
    yp(4) = nx/n ;
    yp(5) = ny/n ;
    yp(6) = nz/n ;
    yp(7) = n;
    yp(8) = sqrt(yp(1)^2+yp(2)^2+yp(3)^2); %path_length
    
    solution = [yp(1); yp(2); yp(3); yp(4); yp(5); yp(6); yp(7); yp(8)];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value,isterminal,direction] = cutoff(t,y,domain)
    
    flag_out=checkifinsidedomain3D([y(1) y(2) y(3)],domain);
    
    if flag_out~=1
        
        value = flag_out;
        isterminal = 1;            % Stop at local minimum
        direction = 0;
        
    else
        
        [n, nx, ny ,nz] = refractive_interpolation(y(1), y(2), y(3),domain);
        
        flag_cutoff=1;
        
        if n<0.001
            flag_cutoff=0;
        end
        
        value = flag_cutoff;
        isterminal = 1;            % Stop at local minimum
        direction = 0;            % [local minimum, local maximum]
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

