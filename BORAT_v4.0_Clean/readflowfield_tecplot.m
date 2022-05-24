function domain=readflowfield_tecplot(flowfieldfilename,f,Cartesian,B1,B2,B3,MagneticField,collisionmodel, shrinkfactor,indexofelectronspecies,Species,pltformat,mixture_name,state_model_str,nb_species,indexof1stspecies,nb_temperatures,indexof1sttemperature,compositiontype,indexoftotaldensity)
    %% this function reads the flowfield and assign the optical properties ri and ac
    
    %%	input
    %%	string			flowfieldname
    %%  double          f                   transmission frequency, Hz
    %%  double          collisionmodel
    %%  double          shrinkfactor
    %%  integer         pltformat
    %%	output
    %%	struct			domain
    %%  internal parameters
    molarmasselectron=5.48579909070e-7;
    nA=6.022140857e23;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% begin input for mutation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    string_number_of_inputs = ['1 2 ', int2str(nb_species),' ', int2str(nb_species+2),' ', int2str(nb_temperatures)];
    prop_mut_names = {'Z'}; % Specific name I gave to the collision frequency in the file mppcalc.cpp, For other specific props check the file
    outcome_length = 1;
    %options.erase_temporary_file = true;
    options.erase_temporary_file = false;
    options.display = true;
    %options.display = false;
    partialdensityfix=1e-30;    %if any partial density is below this value, set it to this value
    temperaturefix=150;         %if any temperature is below this value, set it to this value
    
    %% molar mass of each species
    MM=zeros(nb_species,1);
    
    if  strcmp(mixture_name,'air11Jan')
        % air11
        MM(Species.e)=5.48579909070e-7; MM(Species.N)=7e-3; MM(Species.O)=8e-3; MM(Species.N2)=14e-3; MM(Species.NO)=15e-3; MM(Species.O2)=16e-3; MM(Species.N2p)=14e-3;	MM(Species.NOp)=15e-3;	MM(Species.Np)=14e-3;	MM(Species.O2p)=16e-3;	MM(Species.Op)=8E-3;
        %em                     N           O           N2            NO            O2            N2p               NOp             Np              O2p             Op
    elseif  strcmp(mixture_name,'air11Johannes')
        % air11
        MM(Species.e)=5.48579909070e-7; MM(Species.N)=7e-3; MM(Species.O)=8e-3; MM(Species.N2)=14e-3; MM(Species.NO)=15e-3; MM(Species.O2)=16e-3; MM(Species.N2p)=14e-3;	MM(Species.NOp)=15e-3;	MM(Species.Np)=14e-3;	MM(Species.O2p)=16e-3;	MM(Species.Op)=8E-3;
        %em                             N                   O                   N2                    NO                    O2                    N2p                       NOp                     Np                      O2p                     Op
    elseif   strcmp(mixture_name ,'Mars_CO2_N_ionized')
        % Mars_CO2_N
        MM(Species.e)=5.48579909070e-7; MM(Species.CO2)=44e-3; MM(Species.N2)=28e-3; MM(Species.C)=12e-3; MM(Species.N)=14e-3; MM(Species.O)=16e-3; MM(Species.O2)=32e-3;	MM(Species.CO)=28e-3;	MM(Species.NO)=30e-3;	MM(Species.COp)=28e-3;	MM(Species.NOp)=30E-3;    MM(Species.Cp)=12e-3;	MM(Species.Op)=16e-3;	MM(Species.O2p)=32E-3;
        %e-                     CO2            N2            C            N            O            O2              CO              NO              CO+             NO+               C+            O+              O2+
    elseif   strcmp(mixture_name ,'Mars_CO2_N_ionized_14')
        % Mars_CO2_N_14
        MM(Species.e)=5.48579909070e-7; MM(Species.CO2)=44e-3; MM(Species.N2)=28e-3; MM(Species.C)=12e-3; MM(Species.N)=14e-3; MM(Species.O)=16e-3; MM(Species.O2)=32e-3;	MM(Species.Co)=28e-3;	MM(Species.NO)=30e-3;	MM(Species.COp)=28e-3;	MM(Species.NOp)=30E-3;    MM(Species.Cp)=12e-3;	MM(Species.Op)=16e-3;	MM(Species.O2p)=32E-3;
        %e-                     CO2            N2            C            N            O            O2              CO              NO              CO+             NO+               C+            O+              O2+
    elseif strcmp(mixture_name,'air_11_plato')
        % air11
        MM(Species.e)=5.48579909070e-7; MM(Species.N)=14e-3; MM(Species.O)=16e-3; MM(Species.N2)=28e-3; MM(Species.NO)=30e-3; MM(Species.O2)=32e-3; MM(Species.N2p)=28e-3;	MM(Species.NOp)=30e-3;	MM(Species.Np)=14e-3;	MM(Species.O2p)=32e-3;	MM(Species.Op)=16E-3;
        %em                             N                    O                    N2                    NO                    O2                    N2p                     NOp                     Np                      O2p                     Op
        %       e-  N O N2 NO O2 N2+ NO+ N+ O2+ O+
        
    elseif strcmp(mixture_name,'air_11_ICP')
        % air11
        MM(Species.e)=5.48579909070e-7; MM(Species.N)=14e-3; MM(Species.O)=16e-3; MM(Species.N2)=28e-3; MM(Species.NO)=30e-3; MM(Species.O2)=32e-3; MM(Species.N2p)=28e-3;	MM(Species.NOp)=30e-3;	MM(Species.Np)=14e-3;	MM(Species.O2p)=32e-3;	MM(Species.Op)=16E-3;
        %em                             N                    O                    N2                    NO                    O2                    N2p                     NOp                     Np                      O2p                     Op
        %       e-  N O N2 NO O2 N2+ NO+ N+ O2+ O+
    elseif  strcmp(mixture_name,'7species')
        % air11
        MM(Species.e)=5.48579909070e-7; MM(Species.N)=7e-3; MM(Species.O)=8e-3; MM(Species.N2)=14e-3; MM(Species.NO)=15e-3; MM(Species.O2)=16e-3; MM(Species.N2p)=14e-3;	MM(Species.NOp)=15e-3;
        %em                             N                   O                   N2                    NO                    O2                    N2p                       NOp
    elseif  strcmp(mixture_name,'3species')
        % air11
        MM(Species.e)=5.48579909070e-7; MM(Species.N)=7e-3; MM(Species.O)=8e-3; 
        %em                             N                   O                   
    
    else
        error('mixture not identified')
    end
    %% neutral CO2
    %MM(1)=12e-3;	MM(2)=28e-3;	MM(3)=44e-3;	MM(4)=16e-3;	MM(5)=5.48579909070e-7;
    %%	MM(1)=C;    MM(2)=CO;		MM(3)=CO2;		MM(4)=O;		MM(5)=e-;
    
    
    %    ADD_MIXTURE_PROPERTY(Mm,
    %    "mixture molar mass", "kg/mol", 1,
    %    v[0] = mix().mixtureMw()
    %   )
    %    prop_mut_names2 = {'Mm'};
    %    [MMMM] = retrieve_and_sort([1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1],mixture_name,state_model_str,string_number_of_inputs,prop_mut_names2,outcome_length,options)
    %    size(MMMM)
    %    pause
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% end input for mutation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    z=0;
    id=fopen(flowfieldfilename);
    %read titleline, to be discarded
    fgetl(id);
    %read variables' name line, to be discarded, but number of variables to be
    %determined by number of quotation marks
    
    while ~feof(id)
        z=z+1;
        if pltformat==1
            variablenamesline=fgetl(id);
            domain.nova=size(strfind(variablenamesline,'"'),2)/2;
            dimline=fgetl(id); %% read line that contains grid information
            dimline(strfind(dimline, ' ')) = []            ;
            ipos=strfind(dimline,'i=');jpos=strfind(dimline,'j=');fpos=strfind(dimline,'f=');Npos=strfind(dimline,'N=');Epos=strfind(dimline,'E=');Fpos=strfind(dimline,'F=');
        
        elseif pltformat==11
            variablenamesline=fgetl(id);
            domain.nova=size(strfind(variablenamesline,'"'),2)/2;
            dimline=fgetl(id); %% read line that contains grid information
            dimline(strfind(dimline, ' ')) = []            ;
            ipos=strfind(dimline,'i=');jpos=strfind(dimline,'j=');fpos=strfind(dimline,'f=');Npos=strfind(dimline,'N=');Epos=strfind(dimline,'E=');Fpos=strfind(dimline,'F=');
        
        elseif pltformat==2
            %! read variablesfilenames and count them until key word "ZONE" is read
            domain.nova=0;
            while 1
                if strfind(fgetl(id),'ZONE')
                    break
                end
                domain.nova=domain.nova+1;
            end
            fgetl(id);
            dimline=fgetl(id); %% read line that contains grid information
            dimline(strfind(dimline, ' ')) = []            ;
            ipos=strfind(dimline,'i=');jpos=strfind(dimline,'j=');fpos=strfind(dimline,'f=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos=strfind(dimline,'ZONETYPE=')+7;
            fgetl(id);fgetl(id);
            
        elseif pltformat==3 % This is a STEFANO provided format
            domain.nova=14;
            dimline=fgetl(id); %% read line that contains grid information
            ipos=strfind(dimline,'i=');jpos=strfind(dimline,'j=');fpos=strfind(dimline,'f=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos=strfind(dimline,'ZONETYPE=')+7;
            fgetl(id);fgetl(id);
            
        elseif pltformat==4 % VKI
            %! read variablesfilenames and count them until key word "ZONE" is read
            if z==1
                domain.nova=0;
                while 1
                    if strfind(fgetl(id),'ZONE')
                        break
                    end
                    domain.nova=domain.nova+1;
                end
                fgetl(id);
                dimline=fgetl(id); %% read line that contains grid information
                dimline(strfind(dimline, ' ')) = []            ;
                ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos= [];%Fpos=strfind(dimline,'ZONETYPE=')+7;
                fgetl(id);fgetl(id);
            elseif z>1
                fgetl(id);
                dimline=fgetl(id); %% read line that contains grid information
                dimline(strfind(dimline, ' ')) = []            ;
                ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos= [];%Fpos=strfind(dimline,'ZONETYPE=')+7;
                fgetl(id);fgetl(id);
            else
                error('unknown posformat: %d',pltformat);
            end
            
        elseif pltformat==5 % ICP
            if z==1
                domain.nova=0;
                while 1
                    if strfind(fgetl(id),'ZONE')
                        break
                    end
                    domain.nova=domain.nova+1;
                end
                fgetl(id);
                dimline=fgetl(id); %% read line that contains grid information
                dimline(strfind(dimline, ' ')) = []            ;
                ipos=strfind(dimline,'i=');jpos=strfind(dimline,'j=');fpos=strfind(dimline,'f=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos=strfind(dimline,'ZONETYPE=')+7;
                fgetl(id);fgetl(id);fgetl(id);fgetl(id);fgetl(id);
            else
                z=1;
                break
            end
        elseif pltformat==6 % VKI interpolated
            %! read variablesfilenames and count them until key word "ZONE" is read
            if z==1
                domain.nova=0;
                while 1
                    if strfind(fgetl(id),'ZONE')
                        break
                    end
                    domain.nova=domain.nova+1;
                end
                dimline=fgetl(id); %% read line that contains grid information
                dimline(strfind(dimline, ' ')) = []            ;
                ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos= [];%Fpos=strfind(dimline,'ZONETYPE=')+7;
                fgetl(id);fgetl(id);
            elseif z>1
                fgetl(id);
                dimline=fgetl(id); %% read line that contains grid information
                dimline(strfind(dimline, ' ')) = []            ;
                ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos= [];%Fpos=strfind(dimline,'ZONETYPE=')+7;
                fgetl(id);fgetl(id);
            else
                error('unknown posformat: %d',pltformat);
            end
            
        elseif pltformat==7 % interpolated triangular
            if z==1
            domain.nova=0;
            while 1
                if strfind(fgetl(id),'ZONE')
                    break
                end
                domain.nova=domain.nova+1;
            end
            fgetl(id);
            dimline=fgetl(id); %% read line that contains grid information
            dimline(strfind(dimline, ' ')) = []            ;
            ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos=strfind(dimline,'ZONETYPE=')+7;%Fpos=strfind(dimline,'ZONETYPE=')+7;
            fgetl(id);
            fgetl(id);
            elseif z>1
                fgetl(id);
            dimline=fgetl(id); %% read line that contains grid information
            dimline(strfind(dimline, ' ')) = []            ;
            ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos=strfind(dimline,'ZONETYPE=')+7;%Fpos=strfind(dimline,'ZONETYPE=')+7;
            fgetl(id);
            fgetl(id);
            end
        elseif pltformat==8 % interpolated triangular
            if z==1
            domain.nova=0;
            while 1
                if strfind(fgetl(id),'ZONE')
                    break
                end
                domain.nova=domain.nova+1;
            end
            fgetl(id);
            dimline=fgetl(id); %% read line that contains grid information
            dimline(strfind(dimline, ' ')) = []            ;
            ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos=strfind(dimline,'ZONETYPE=')+7;%Fpos=strfind(dimline,'ZONETYPE=')+7;
            fgetl(id);
            fgetl(id);
            fgetl(id);
            fgetl(id);
           
         
            elseif z>1
                fgetl(id);
            dimline=fgetl(id); %% read line that contains grid information
            dimline(strfind(dimline, ' ')) = []            ;
            ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos=strfind(dimline,'ZONETYPE=')+7;%Fpos=strfind(dimline,'ZONETYPE=')+7;
            fgetl(id);
            fgetl(id);
            fgetl(id);
            fgetl(id);
            end
      elseif pltformat==9 % VKI
            %! read variablesfilenames and count them until key word "ZONE" is read
            if z==1
                domain.nova=0;
                while 1
                    if strfind(fgetl(id),'ZONE')
                        break
                    end
                    domain.nova=domain.nova+1;
                end
                fgetl(id);
                dimline=fgetl(id); %% read line that contains grid information
                dimline(strfind(dimline, ' ')) = []            ;
                ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos=strfind(dimline,'ZONETYPE=')+7;
                fgetl(id);fgetl(id);
            elseif z>1
                fgetl(id);
                dimline=fgetl(id); %% read line that contains grid information
                dimline(strfind(dimline, ' ')) = []            ;
                ipos=strfind(dimline,'I=');jpos=strfind(dimline,'J=');fpos=strfind(dimline,'K=');Npos=strfind(dimline,'Nodes=')+4;Epos=strfind(dimline,'Elements=')+7;Fpos=strfind(dimline,'ZONETYPE=')+7;
                fgetl(id);fgetl(id);
            end
        else
            error('unknown pltformat: %d',pltformat);
        end
        
        if not(isempty(ipos)) && not(isempty(jpos)) && not(isempty(fpos)) && isempty(Npos) && isempty(Epos) && isempty(Fpos)
            domain.unstructured=0;
            domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(ipos+2:jpos-1))*str2double(dimline(jpos+2:fpos-1));
            domain.( strcat('zone',num2str(z)) ).E=0;
            fprintf('\t\tstructured mesh, zone=%d: i=%d j=%d \n',z,str2double(dimline(ipos+2:jpos-1)),str2double(dimline(jpos+2:fpos-1)));
        elseif isempty(ipos) && isempty(jpos) && isempty(fpos) && not(isempty(Npos)) && not(isempty(Epos)) && not(isempty(Fpos))
            domain.unstructured=1;
            if pltformat==1
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2));
            elseif pltformat==2
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1-7));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2-7));
            elseif pltformat==3
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1-7));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2-7));
            elseif pltformat==4
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1-7));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2-7));
            elseif pltformat==5
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1-7));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2-7));
            elseif pltformat==6
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1-7));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2-7));
            elseif pltformat==7
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1-7));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2-7));
            elseif pltformat==8
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1-7));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2-7));
            elseif pltformat==9
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1-7));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2-7));
            elseif pltformat==11
                domain.( strcat('zone',num2str(z)) ).N=str2double(dimline(Npos+2:Epos-1));
                domain.( strcat('zone',num2str(z)) ).E=str2double(dimline(Epos+2:Fpos-2));
            end
        else
            fprintf('\n ipos %d jpos %s fpos %s Npos %s Epos %s Fpos %s',ipos,jpos,fpos,Npos,Epos,Fpos);
            error('neither structured nor unstructured mesh recognised. \n FYI: pltformat=%d',pltformat)
        end
        %%read data, call function to compute optical properties
        domain.( strcat('zone',num2str(z)) ).variables=zeros(  domain.nova  ,  domain.( strcat('zone',num2str(z)) ).N);
        for i=1:domain.( strcat('zone',num2str(z)) ).N
            line=fgetl(id);
            domain.( strcat('zone',num2str(z)) ).variables(1:domain.nova,i)=sscanf(line,'  %g');
            
            %%compute numberdensity
            if compositiontype==1
                %  domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i)*domain.( strcat('zone',num2str(z)) ).variables(indexoftotaldensity,i)*nA/molarmasselectron;
                % the formula has been corrected, error in composition type
                domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i)*nA/molarmasselectron;
                
            elseif compositiontype==2
                domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i);
                
            elseif compositiontype==3
                MassFrac= zeros(nb_species,1);
                for Spec = 1: nb_species
                    MassFrac(Spec) = (domain.( strcat('zone',num2str(z)) ).variables((indexof1stspecies+Spec-1),i)* MM(Spec));
                end
                MassFrac_tot = sum(MassFrac);
                domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexofelectronspecies,i)*domain.( strcat('zone',num2str(z)) ).variables(indexoftotaldensity,i)*nA/MassFrac_tot;
                %the formula has been corrected, error in composition type
                %domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i)*nA*/molarmasselectron;
                
            else
                error('unknown compositiontype: %d',compositiontype);
            end
            
            %%collisionfrequency placeholder
            domain.( strcat('zone',num2str(z)) ).variables(domain.nova+2,i)  = 0;
        end
        if pltformat ==11 % removes i variable from IRS solution
            domain.( strcat('zone',num2str(z)) ).variables=domain.( strcat('zone',num2str(z)) ).variables(2:end,:);
            domain.nova=domain.nova-1;
        end
        Q=domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,:);
                Q(isnan(Q))=0;
                domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,:)=Q;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%begin mutation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if collisionmodel==1
            %% assignment of species for mutation
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
            
            %% assignment of temperatures for mutation
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
            %% call mutation interface
            fprintf('\n \nmutation++ begin \n');
            [collfreq] = retrieve_and_sort(fields,mixture_name,state_model_str,string_number_of_inputs,prop_mut_names,outcome_length,options);
            
            fprintf('mutation++ end \n \n');
            %% compute collision frequency
            domain.( strcat('zone',num2str(z)) ).variables(domain.nova+2,:)  = collfreq';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%end mutation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %     domain.nova+1 -> electron density
        %     domain.nova+2 -> coll freq
        %     domain.nova+3 -> X
        %     domain.nova+4 -> Z
        %     domain.nova+5 -> ri
        %     domain.nova+6 -> ac
       
        
        % needed to initialize the domain in case of composition equal to
        % zero
        for i=1:domain.( strcat('zone',num2str(z)) ).N
            
            domain.(strcat('zone',num2str(z))).variables(domain.nova+3,i)=0;    %X
            domain.(strcat('zone',num2str(z))).variables(domain.nova+4,i) =0;   %Z
            domain.(strcat('zone',num2str(z))).variables(domain.nova+5,i)=0;    %ri
            domain.(strcat('zone',num2str(z))).variables(domain.nova+6,i)=0;    %ac
            AppletonMagnetY(1,i)=0;
            AppletonMagnetB(1,i)=0; 
        end
        %
        Nova=domain.( strcat('nova'));
        for i=1:domain.( strcat('zone',num2str(z)) ).N
            %% optical properties, i.e. refractive index and absorption coefficient
            [domain.(strcat('zone',num2str(z))).variables(domain.nova+3,i),domain.(strcat('zone',num2str(z))).variables(domain.nova+4,i),domain.(strcat('zone',num2str(z))).variables(domain.nova+5,i),domain.(strcat('zone',num2str(z))).variables(domain.nova+6,i)]  =  opticalpropertiesMagnetic(f,collisionmodel,Cartesian,B1,B2,B3,MagneticField, domain.(strcat('zone',num2str(z))).variables(:,i),Nova );
        end
        %from here absorption is in db
        pause
        %% read elements (this necessary to proceed to the next zone)
        
        
        if pltformat==7||pltformat==8||pltformat==11 %triangular mesh
            b_point=0;
            domain.( strcat('zone',num2str(z)) ).elements=zeros(domain.( strcat('zone',num2str(z)) ).E,3 );
        else
            domain.( strcat('zone',num2str(z)) ).elements=zeros(domain.( strcat('zone',num2str(z)) ).E,4 );
        end
        
        for i=1:domain.( strcat('zone',num2str(z)) ).E
            
            line=fgetl(id);
            domain.( strcat('zone',num2str(z)) ).elements(i,:)=sscanf(line,'  %d');
            
            
            if pltformat==7||pltformat==8||pltformat==11
                
                if domain.( strcat('zone',num2str(z)) ).elements(i,1)==domain.( strcat('zone',num2str(z)) ).elements(i,2)
                    b_point=b_point+1;
                end
            end
            
        end
        
        if pltformat==7||pltformat==8||pltformat==11
            
            domain.( strcat('zone',num2str(z)) ).NodeMap=zeros(domain.( strcat('zone',num2str(z)) ).E-b_point,3);
            domain.( strcat('zone',num2str(z)) ).Boundary_points=zeros(b_point,2);
            
            
            for i=1:domain.( strcat('zone',num2str(z)) ).E
                if domain.( strcat('zone',num2str(z)) ).elements(i,1)==domain.( strcat('zone',num2str(z)) ).elements(i,2)
                    domain.( strcat('zone',num2str(z)) ).Boundary_points(i,1)= domain.( strcat('zone',num2str(z)) ).elements(i,1);
                    domain.( strcat('zone',num2str(z)) ).Boundary_points(i,2)= domain.( strcat('zone',num2str(z)) ).elements(i,2);
                else
                    domain.( strcat('zone',num2str(z)) ).NodeMap(i-b_point,:)=domain.( strcat('zone',num2str(z)) ).elements(i,:);
                end
            end
            
        end
        
        if pltformat==7||pltformat==8||pltformat==11 % for TRI mesh
            
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
            
        else %quadmesh
    %        
    %        % define boundary and volume of zone
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
            %% New Meshing thechnique
        
      %  domain.( strcat('zone',num2str(z)) ).shp=alphaShape(domain.( strcat('zone',num2str(z)) ).variables(1,:)',domain.( strcat('zone',num2str(z)) ).variables(2,:)');
      %  a = criticalAlpha(domain.( strcat('zone',num2str(z)) ).shp,'all-points')*1.05;
      %  domain.( strcat('zone',num2str(z)) ).shp=alphaShape(domain.( strcat('zone',num2str(z)) ).variables(1,:)',domain.( strcat('zone',num2str(z)) ).variables(2,:)',a);
      %  domain.( strcat('zone',num2str(z)) ).shp2=alphaTriangulation(domain.( strcat('zone',num2str(z)) ).shp);
       % domain.( strcat('zone',num2str(z)) ).delaunay=domain.( strcat('zone',num2str(z)) ).shp2;
        end
        dimline=fgetl(id);

        
      
        %% Boundary from Geometry from Mesh
        % Create geometry
        nodes = [domain.( strcat('zone',num2str(z)) ).variables(1,:);domain.( strcat('zone',num2str(z)) ).variables(2,:)];
        elements = domain.( strcat('zone',num2str(z)) ).delaunay';
        model = createpde();
        [domain.( strcat('zone',num2str(z)) ).Geometry,domain.( strcat('zone',num2str(z)) ).meshSolution]=geometryFromMesh(model,nodes,elements);
        %Find nodes on edges
        for nodals = 1 : domain.( strcat('zone',num2str(z)) ).Geometry.NumEdges
            domain.( strcat('zone',num2str(z)) ).nf.( strcat('N',num2str(nodals)) ) = findNodes(domain.( strcat('zone',num2str(z)) ).meshSolution,'region','Edge',nodals);
        end
        %Save node vectors in cell
        domain.( strcat('zone',num2str(z)) ).nf.cell=struct2cell(domain.( strcat('zone',num2str(z)) ).nf);
        %Get first and last edge vctor element
        domain.( strcat('zone',num2str(z)) ).nf.boundary=[ cellfun(@(v)v(1),domain.( strcat('zone',num2str(z)) ).nf.cell),cellfun(@(v)v(end),domain.( strcat('zone',num2str(z)) ).nf.cell)];
        %Sort edge vectors
        SortVector=[];
        SortVector(1)=1;
        SortVector(2)=find(domain.( strcat('zone',num2str(z)) ).nf.boundary(:,1)==domain.( strcat('zone',num2str(z)) ).nf.boundary(1,2));
        for i=3:domain.( strcat('zone',num2str(z)) ).Geometry.NumEdges
        SortVector(i)=find(domain.( strcat('zone',num2str(z)) ).nf.boundary(:,1)==domain.( strcat('zone',num2str(z)) ).nf.boundary(SortVector(i-1),2));
        end
        domain.( strcat('zone',num2str(z)) ).nf.cellsorted=domain.( strcat('zone',num2str(z)) ).nf.cell(SortVector);
        %Combine al edge vectors
        domain.( strcat('zone',num2str(z)) ).nf.final=[domain.( strcat('zone',num2str(z)) ).nf.cellsorted{:}];
        %Save nodes in bound
        domain.( strcat('zone',num2str(z)) ).bound=domain.( strcat('zone',num2str(z)) ).nf.final;
        %Plot geometry and boundary
       % figure
     %   pdegplot(model,'FaceLabels','on','EdgeLabels','on','FaceAlpha',0.5)
        
     %   figure
     %   plot (domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:)) )
        
        %end boundary from geometry from mesh
       
    end % while eof
    fclose(id);
    domain.nozones=z;
    domain.nova=domain.nova+6; %%make this automatic ???
end %% function read_flowfield