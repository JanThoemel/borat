function domain=import_flowfield_tecplot(flowfieldfilename,f,collisionmodel, shrinkfactor,pltformat,mixture_name,state_model_str,nb_species,indexof1stspecies,nb_temperatures,indexof1sttemperature,compositiontype,indexoftotaldensity)
    %%% this function import the szplt file with the flowfield and assign the optical properties ri and ac
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  TecIO Setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    tecplot_home = '/Applications/Tecplot 360 EX 2020 R2/Tecplot 360 EX 2020 R2.app/Contents';
    tecio_path = strcat(tecplot_home, '/Frameworks/libtecio.dylib');
    % tecio_header_path = '/Applications/Tecplot 360 EX 2020 R2/Tecplot 360 EX 2020 R2.app/Contents/include/TECIO.h'; %Uses simplified TECIO.h included in package.
    tecio_header_path = '/TecIO/TECIO.h'; %Uses simplified TECIO.h included in package.
    
    
    % tecplot_home = '/apps/leuven/skylake/2018a/software/Tecplot/2018r2/360ex_2018r2';
    % tecio_path = strcat(tecplot_home, '/bin/libtecio.so');
    % tecio_header_path = '/TecIO/TECIO.h'; %Uses simplified TECIO.h included in package.
    
    if ~libisloaded('tecio')
        [notfound, warnings] = loadlibrary(tecio_path, tecio_header_path,...
            'alias', 'tecio');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% begin input for mutation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %  internal parameters
    molarmasselectron=5.48579909070e-7;
    nA=6.022140857e23;
    
    string_number_of_inputs = ['1 2 ', int2str(nb_species),' ', int2str(nb_species+2),' ', int2str(nb_temperatures)];
    prop_mut_names = {'Z'}; % Specific name I gave to the collision frequency in the file mppcalc.cpp, For other specific props check the file
    outcome_length = 1;
    
    options.erase_temporary_file = false;
    options.display = true;
    
    partialdensityfix=1e-30;    %if any partial density is below this value, set it to this value
    temperaturefix=150;         %if any temperature is below this value, set it to this value
    
    % molar mass of each species
    MM=zeros(nb_species,1);
    if  strcmp(mixture_name,'air_11')
        % air11
        MM(1)=5.48579909070e-7; MM(2)=7e-3; MM(3)=8e-3; MM(4)=14e-3; MM(5)=15e-3; MM(6)=16e-3; MM(7)=14e-3;	MM(8)=15e-3;	MM(9)=14e-3;	MM(10)=16e-3;	MM(11)=8E-3;
        %em                      N           O           N2           NO           O2           N2p              NOp             Np              O2p             Op
    elseif strcmp(mixture_name,'air_11_plato')
        % air11
        MM(1)=5.48579909070e-7; MM(2)=14e-3; MM(3)=16e-3; MM(4)=28e-3; MM(5)=30e-3; MM(6)=32e-3; MM(7)=28e-3;	MM(8)=30e-3;	MM(9)=14e-3;	MM(10)=32e-3;	MM(11)=16E-3;
        %em                      N           O           N2           NO           O2           N2p              NOp             Np              O2p             Op
    elseif strcmp(mixture_name,'air_11_ICP')
        % air11
        MM(1)=5.48579909070e-7; MM(2)=14e-3; MM(3)=16e-3; MM(4)=28e-3; MM(5)=30e-3; MM(6)=32e-3; MM(7)=28e-3;	MM(8)=30e-3;	MM(9)=14e-3;	MM(10)=32e-3;	MM(11)=16E-3;
        %em                      N           O           N2           NO           O2           N2p              NOp             Np              O2p             Op
    elseif   strcmp(mixture_name ,'Mars_CO2_N_ionized')
        % Mars_CO2_N
        MM(1)=5.48579909070e-7; MM(2)=44e-3; MM(3)=28e-3; MM(4)=12e-3; MM(5)=14e-3; MM(6)=16e-3; MM(7)=32e-3;	MM(8)=28e-3;	MM(9)=30e-3;	MM(10)=28e-3;	MM(11)=30E-3;    MM(12)=12e-3;	MM(13)=16e-3;	MM(14)=32E-3;
        %e-                      CO2         N2           C            N            O            O2             CO              NO              CO+             NO+             C+              O+              O2+
    elseif   strcmp(mixture_name ,'Mars_CO2_N_ionized_14')
        % Mars_CO2_N_14
        MM(1)=5.48579909070e-7; MM(2)=44e-3; MM(3)=28e-3; MM(4)=12e-3; MM(5)=14e-3; MM(6)=16e-3; MM(7)=32e-3;	MM(8)=28e-3;	MM(9)=30e-3;	MM(10)=28e-3;	MM(11)=30E-3;    MM(12)=12e-3;	MM(13)=16e-3;	MM(14)=32E-3;
        %e-                      CO2         N2           C            N            O            O2             CO              NO              CO+             NO+             C+              O+              O2+
        
    else
        error('mixture not identified')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % end input for mutation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Open szplt file
    
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
        %% Get information on what the Zone Type is.
        
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
            
            [isok, ~, domain.( strcat('zone',num2str(z)) ).variables(var,:)] = calllib('tecio', 'tecZoneVarGetFloatValues', ...
                handle, z,var, 1,numValues,domain.( strcat('zone',num2str(z)) ).variables(var,:));
            
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
                
                
                if domain.( strcat('zone',num2str(z)) ).zoneType==2 % for TRI mesh
                    
                    b_point=0;
                    for i=1:3:numValues
                        
                        if nodeMap(i)==nodeMap(i+1)
                            b_point=b_point+1;
                        end
                        
                    end
                    
                    domain.( strcat('zone',num2str(z)) ).NodeMap=zeros(jMax-b_point,3);
                    domain.( strcat('zone',num2str(z)) ).Boundary_szplt=zeros(b_point,2);
                    
                    k=1;
                    b=1;
                    for i=1:3:numValues
                        if nodeMap(i)==nodeMap(i+1)
                            domain.( strcat('zone',num2str(z)) ).Boundary_szplt(b,1)= nodeMap(i+1);
                            domain.( strcat('zone',num2str(z)) ).Boundary_szplt(b,2)= nodeMap(i+2);
                            b=b+1;
                        else
                            domain.( strcat('zone',num2str(z)) ).NodeMap(k,1)=nodeMap(i);
                            domain.( strcat('zone',num2str(z)) ).NodeMap(k,2)=nodeMap(i+1);
                            domain.( strcat('zone',num2str(z)) ).NodeMap(k,3)=nodeMap(i+2);
                            k=k+1;
                        end
                    end
                    
                end
            end
            
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% begin mutation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i=1:domain.( strcat('zone',num2str(z)) ).N
            
            %%compute numberdensity
            if compositiontype==1 %check this
                %  domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i)*domain.( strcat('zone',num2str(z)) ).variables(indexoftotaldensity,i)*nA/molarmasselectron;
                domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i)*nA/molarmasselectron;
                
            elseif compositiontype==2
                domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i);
                
            elseif compositiontype==3
                
                Mtot=0;
                for k=1:nb_species
                    Mtot=Mtot+MM(k)*domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies-1+k,i);
                end
                
                domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i)*domain.( strcat('zone',num2str(z)) ).variables(indexoftotaldensity,i)*nA...
                    /(Mtot);
                
            else
                error('unknown compositiontype: %d',compositiontype);
            end
            %%collisionfrequency placeholder
            domain.( strcat('zone',num2str(z)) ).variables(domain.nova+2,i)  = 0;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if collisionmodel==1
            %%assignment of species for mutation
            fields(:,1:nb_species) = domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies:indexof1stspecies+nb_species-1,:)';
            if compositiontype==2 %% if data is available as number density [1/m3], convert to partial densities [kg/m3]
                for i=1:size(fields,1)
                    for j=1:size(fields,2)
                        fields(i,j)=fields(i,j)/nA*MM(j);
                    end
                end
            end
            
            if compositiontype==3 %% if data is available as mass fraction, convert to partial densities [kg/m3]
                for i=1:size(fields,1)
                    
                    Mtot=0;
                    for j=1:size(fields,2)
                        Mtot=Mtot+MM(j)*fields(i,j);
                    end
                    
                    for j=1:size(fields,2)
                        fields(i,j)=fields(i,j)*MM(j)/Mtot*domain.( strcat('zone',num2str(z)) ).variables(indexoftotaldensity,i);
                    end
                    
                end
            end
            
            for i=1:size(fields,1) %%check if partial densities might be to low, if so fix it, dont fix electron density
                for j=1:size(fields,2)
                    if fields(i,j)<partialdensityfix
                        fields(i,j)=partialdensityfix;
                    end
                end
            end
            
            % assignment of temperatures for mutation
            fields(:,nb_species+1:nb_species+1-1+nb_temperatures) = domain.( strcat('zone',num2str(z)) ).variables(indexof1sttemperature:indexof1sttemperature-1+nb_temperatures,:)';
            for i=1:size(fields,1)
                %fprintf('\n m1%f',domain.( strcat('zone',num2str(z)) ).variables(indexof1sttemperature-1:indexof1sttemperature-1+nb_temperatures-1,i));
                %fprintf('\n   %f',domain.( strcat('zone',num2str(z)) ).variables(indexof1sttemperature:indexof1sttemperature-1+nb_temperatures,i));
                %fprintf('\n p1%f',domain.( strcat('zone',num2str(z)) ).variables(indexof1sttemperature+1:indexof1sttemperature-1+nb_temperatures+1,i));
                %pause
                if fields(i,nb_species+1)<temperaturefix
                    fields(i,nb_species+1)=temperaturefix;
                end
            end
            %if (nb_species+3)*domain.( strcat('zone',num2str(z)) ).N > 1e6
            %    fprintf('\n  nb_species  %d  number of nodes  %d  (nb_species+3)*number of nodes %d  \n', nb_species,domain.( strcat('zone',num2str(z)) ).N,(nb_species+3)*domain.( strcat('zone',num2str(z)) ).N);
            %    pause
            %end
            % call mutation interface
            fprintf('\n \nmutation++ begin \n');
            [collfreq] = retrieve_and_sort(fields,mixture_name,state_model_str,string_number_of_inputs,prop_mut_names,outcome_length,options);
            
            fprintf('mutation++ end \n \n');
            %compute collision frequency
            domain.( strcat('zone',num2str(z)) ).variables(domain.nova+2,:)  = collfreq';
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%end mutation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% compute optical properties
        
        %     domain.nova+1 -> electron density
        %     domain.nova+2 -> coll freq
        
        for i=1:domain.( strcat('zone',num2str(z)) ).N
            
            domain.(strcat('zone',num2str(z))).variables(domain.nova+3,i)=0;    %X
            domain.(strcat('zone',num2str(z))).variables(domain.nova+4,i) =0;   %Z
            domain.(strcat('zone',num2str(z))).variables(domain.nova+5,i)=0;    %ri
            domain.(strcat('zone',num2str(z))).variables(domain.nova+6,i)=0;    %ac
            
        end
        
        
        
        for i=1:domain.( strcat('zone',num2str(z)) ).N
            % optical properties, i.e. refractive index and absorption coefficient
            [domain.(strcat('zone',num2str(z))).variables(domain.nova+3,i),domain.(strcat('zone',num2str(z))).variables(domain.nova+4,i),domain.(strcat('zone',num2str(z))).variables(domain.nova+5,i),domain.(strcat('zone',num2str(z))).variables(domain.nova+6,i)]  =  opticalproperties(f,collisionmodel, domain.(strcat('zone',num2str(z))).variables(:,i) );
        end
        
        %from here absorption is in db
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        pause
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% delauney
        
        
        if domain.( strcat('zone',num2str(z)) ).zoneType==2 % for TRI mesh
            
            domain.( strcat('zone',num2str(z)) ).DT_triangulation=triangulation(domain.( strcat('zone',num2str(z)) ).NodeMap,domain.( strcat('zone',num2str(z)) ).variables(1:2,:)');
            
            domain.( strcat('zone',num2str(z)) ).C=domain.( strcat('zone',num2str(z)) ).DT_triangulation.freeBoundary;
            
            domain.( strcat('zone',num2str(z)) ).bound=zeros(size(domain.( strcat('zone',num2str(z)) ).C,1)+1,1);
            
            for i=1:size(domain.( strcat('zone',num2str(z)) ).bound)-1
                domain.( strcat('zone',num2str(z)) ).bound(i)=domain.( strcat('zone',num2str(z)) ).C(i,1);
            end
            domain.( strcat('zone',num2str(z)) ).bound(end)=domain.( strcat('zone',num2str(z)) ).C(end,2);
            
            
            [~, domain.( strcat('zone',num2str(z)) ).v]=boundary(domain.( strcat('zone',num2str(z)) ).variables(1:2,domain.( strcat('zone',num2str(z)) ).bound(:))',shrinkfactor);
            
            domain.( strcat('zone',num2str(z)) ).DT=delaunayTriangulation(domain.( strcat('zone',num2str(z)) ).variables(1:2,:)',domain.( strcat('zone',num2str(z)) ).C);
            
            tf=isInterior(domain.( strcat('zone',num2str(z)) ).DT);
            % remove all exterior triangles
            counter=1;
            for i=1:size(domain.( strcat('zone',num2str(z)) ).DT,1)
                if tf(i)
                    domain.( strcat('zone',num2str(z)) ).delaunay(counter,:)=domain.( strcat('zone',num2str(z)) ).DT(i,:);
                    counter=counter+1;
                end
            end
            
        elseif domain.( strcat('zone',num2str(z)) ).zoneType==3 % for QUAD mesh
            
            %% define boundary and volume of zone
            [domain.( strcat('zone',num2str(z)) ).bound , domain.( strcat('zone',num2str(z)) ).v]=boundary(domain.( strcat('zone',num2str(z)) ).variables(1:2,:)',shrinkfactor);
            %%define constraints for delaunay triangulation
            domain.( strcat('zone',num2str(z)) ).C(1,:) = [domain.( strcat('zone',num2str(z)) ).bound(size(domain.( strcat('zone',num2str(z)) ).bound,1)) , domain.( strcat('zone',num2str(z)) ).bound(1)];
            for i=2:size(domain.( strcat('zone',num2str(z)) ).bound)
                domain.( strcat('zone',num2str(z)) ).C(i,:)=[domain.( strcat('zone',num2str(z)) ).bound(i-1),domain.( strcat('zone',num2str(z)) ).bound(i)];
            end
            %% define constraint delaunay triangulation
            domain.( strcat('zone',num2str(z)) ).DT=delaunayTriangulation(domain.( strcat('zone',num2str(z)) ).variables(1:2,:)',domain.( strcat('zone',num2str(z)) ).C);
            tf=isInterior(domain.( strcat('zone',num2str(z)) ).DT);
            %% remove all exterior triangles
            counter=1;
            for i=1:size(domain.( strcat('zone',num2str(z)) ).DT,1)
                if tf(i)
                    domain.( strcat('zone',num2str(z)) ).delaunay(counter,:)=domain.( strcat('zone',num2str(z)) ).DT(i,:);
                    counter=counter+1;
                end
            end
            
            
        else
            error('unknown zoneType: %d',domain.( strcat('zone',num2str(z)) ).zoneType);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end %cicle in zone
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    domain.nova=domain.nova+6; %%make this automatic ???
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    calllib('tecio', 'tecFileReaderClose', handle) ;
    
    if libisloaded('tecio')
        unloadlibrary('tecio')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end %% function read_flowfield


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%------