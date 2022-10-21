function        export_flowfield_tecplot(domain,output_filename,collisionmodel,criticaldensity,max_edensity)
    
    
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
    
    if collisionmodel==1
        domain.VarNamesList=strcat(domain.VarNamesList,',','Xe',',','Coll_freq',',','X',',','Z',',','Ri',',','k');
    else
        domain.VarNamesList=strcat(domain.VarNamesList,',','Xe',',','Coll_freq',',','X',',','Z',',','Ri',',','k');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fileFormat=1; %szplt
    data_type=0; %double, gets overwritten after
    
    
    [isok,~,~,~,~,handle_output] = calllib('tecio', 'tecFileWriterOpen', ...
        output_filename, domain.dataTitle, domain.VarNamesList, ...
        fileFormat, domain.fileType, data_type, [],[] );
    
    if isok~=0
        disp('error');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for z=1:domain.nozones
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        for i=(size(domain.( strcat('zone',num2str(z)) ).varTypes,1)+1):domain.nova
            domain.( strcat('zone',num2str(z)) ).varTypes(i,1)=domain.( strcat('zone',num2str(z)) ).varTypes(1,1);
            domain.( strcat('zone',num2str(z)) ).passiveVarList(i,1)=domain.( strcat('zone',num2str(z)) ).passiveVarList(1,1);
            domain.( strcat('zone',num2str(z)) ).valueLocation(i,1)=domain.( strcat('zone',num2str(z)) ).valueLocation(1,1);
            domain.( strcat('zone',num2str(z)) ).shareVarFromZone(i,1)=domain.( strcat('zone',num2str(z)) ).shareVarFromZone(1,1);
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        outputZone = 0;
        
        
        [isok, ~, ~, ~, ~ , ~ , ~ , outputZone] = calllib('tecio', 'tecZoneCreateFE',...
            handle_output, domain.( strcat('zone',num2str(z)) ).zoneTitle, ...
            domain.( strcat('zone',num2str(z)) ).zoneType, ...
            domain.( strcat('zone',num2str(z)) ).iMax, ...
            domain.( strcat('zone',num2str(z)) ).jMax, ...
            domain.( strcat('zone',num2str(z)) ).varTypes, ...
            domain.( strcat('zone',num2str(z)) ).shareVarFromZone, ...
            domain.( strcat('zone',num2str(z)) ).valueLocation, ...
            domain.( strcat('zone',num2str(z)) ).passiveVarList, ...
            domain.( strcat('zone',num2str(z)) ).shareConnectivityFromZone, ...
            domain.( strcat('zone',num2str(z)) ).numFaceConnections, ...
            domain.( strcat('zone',num2str(z)) ).faceNeighborMode,outputZone );
        
        if isok~=0
            disp('error');
        end
        
        
        for var = 1:domain.nova
            
            [isok, ~, domain.( strcat('zone',num2str(z)) ).variables(var,:)] = calllib('tecio', 'tecZoneVarWriteFloatValues',...
                handle_output, outputZone, var, 0, domain.( strcat('zone',num2str(z)) ).N, domain.( strcat('zone',num2str(z)) ).variables(var,:) );
            
            if isok~=0
                disp('error');
            end
            
        end
        
        
        if (domain.( strcat('zone',num2str(z)) ).zoneType ~= 0 && domain.( strcat('zone',num2str(z)) ).shareConnectivityFromZone == 0)
            
            
            if domain.( strcat('zone',num2str(z)) ).is64Bit
                
                
                [isok, ~, domain.( strcat('zone',num2str(z)) ).originalNodeMap(:,1)] = calllib('tecio', 'tecZoneNodeMapWrite64',...
                    handle_output, outputZone,0,1, domain.( strcat('zone',num2str(z)) ).E, domain.( strcat('zone',num2str(z)) ).originalNodeMap(:,1));
                
                if isok~=0
                    disp('error');
                end
                
            else
                
                [isok, ~, domain.( strcat('zone',num2str(z)) ).originalNodeMap(:,1)] = calllib('tecio', 'tecZoneNodeMapWrite32',...
                    handle_output, outputZone,0,1, domain.( strcat('zone',num2str(z)) ).E, domain.( strcat('zone',num2str(z)) ).originalNodeMap(:,1));
                if isok~=0
                    disp('error');
                end
                
            end
            
            
        end
        
    end
    
    
    calllib('tecio', 'tecFileWriterClose', handle_output) ;
    
    if libisloaded('tecio')
        unloadlibrary('tecio')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %end function
