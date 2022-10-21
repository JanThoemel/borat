clear all; clc;format long;close all; initime=cputime;
pause('off'); %there are pause points in the script, enable this option and they won't be active
%pause('on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NAME='BlackOut RAy-Tracer';NAMESHORT='BORAT ';VERSION=2.0;


%% Version 2.0 - March 2022

%% 0. TestCase Definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Raytracing solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solverEikonal=0;
solverSnell=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialisation Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Read Flowfield Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    flowfieldfilename = 'C:/folder/*.dat'; %file path
  
%%%% Define File Format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ASCII or SZPLT File Format (0 inactive - 1 active) %%%%%%%%%%%%%%%%%%%

    read = 1;   % if 1: reads ASCII file (import must be 0 then)
    import = 0; % if 1: imports szplt file (read must be 0 then)


%%%% Export Domain as SZPLT file (Tecplot format) %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    flag_exportdomain = 0; % if 1: export; if 0: no export)
    output_filename = 'flowfield_output.szplt';
     

%%%% Choose ASCII Tecplot File Format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Provided by:                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pltformat= 1; % 1 to 11 dependent on your input file
      

%%%% Create Folder for Run Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    runFolder='RunOutput';
    %runFolder='Debug';

    mkdir (runFolder);

%%%% Save Simulation Solution Files in .mat Format %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % solutionFile='raytracing_solution';
    solutionFile='raytracing_solution';

%%%% Select Transmission Trequency [Hz] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f = 400.0e6;                  %%: UHF 
    
    

%%%% Axis-Symmetry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    symmetryline=1; %%: if 1: symmetry line at y=0; if 0: no symmetry or symmetryline not at y=0 ????


%%%% Choose Shrinkfactor for the matlab boundary function %%%%%%%%%%%%%%%%%
%%%% to determine concave hull of each zone               %%%%%%%%%%%%%%%%%

    %shrinkfactor=0.999;          %%: 
    

%%%% wall boundary limits box (defines wall reflections zone) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case
     limitx1=0; % lower x limit
     limitx2=1; % upper x limit
     
     limity1=-1; % lower y limit
     limity2=1; % upper y limit

 
    
 
%%%% Rays starting points (Antenna Position) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    pooo1=[0.1 0.01];                     %%: 
    pooo2=[0.1 0.01];                     %%: 
    
    
%%%% Min and Max Angle of Antenna Radiation Field [degrees] %%%%%%%%%%%%%%%
%%%% dir1 shall be counterclockwise to dir2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    dir1=90+80;					%%: Max initial radiation angle
    dir2=90-80;					%%: Min initial radiation angle


%%%% Plot LOS to Relay Satellite %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plotscdir=0; %%: if 1: LOS will be plotted


%%%% Spacecraft Direction/Cone Angle [degree] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%% Case
    scdir=90;                  %%:  Direction to s/c
    

%%%% spacecraft direction range(angular) [degree] %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    scdir_range=30;


%%%% Number of Rays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% choose value that allows integer step size, i.e. 179,89,59 etc %%%%%%%

    maxangles=50; %number of rays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Stepsize [m]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stepsize = 10 * CFD cell size (at relevant locations) %%%%%%%%%%%%%%%%

    ss=0.002;   %stepsize for marching in space method


%%%% Max Number of Steps for Ray-Tracing to Escape the Plasma %%%%%%%%%%%%%
%%%% ss*maxsteps > 2*size of domain (rays reach boundary) %%%%%%%%%%%%%%%%%

    maxsteps=5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Physics Modeling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Collsion Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    collisionmodel=0;       %%: (0: no collisions, 1: mutation (Linux))


%%%% Absorptionlimits for Ploting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Absorptionlimits(2)-10 is used to interrupt raytracing because it is %
%%%% assumed the signal strength is too low to continue %%%%%%%%%%%%%%%%%%%

    absorptionlimits=[1e-6 40];


%%%% Number of Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %   nb_species=10;      %%: Number of species


    
%%%% Order of Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Order = zeros(nb_species,1);

%%%% Mars ionized

%Species.e=1; Species.CO2=2; Species.N2=3; Species.C=4; Species.N=5; Species.O=6; Species.O2=7; Species.CO=8; Species.NO=9; Species.COp=10; Species.NOp=11; Species.Cp=12; Species.Op=13; Species.O2p=14;


%%%% Index of First Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    indexof1stspecies=3;    %%: Index of 1st species in input file

    
%%%% Index of Electron Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    indexofelectronspecies=3;    %%: Index of electron species in input file


%%%% Number of Temperatures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nb_temperatures=1;      %%: CNEQ 
    %nb_temperatures=2;      %%: TCNEQ 


%%%% Index of First Temperature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    indexof1sttemperature=19;   %%: Index of 1st temperature in input file

  

%%%% Index of Total Density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    indexoftotaldensity=20;     %%: Index of total density in input file

    
    
%%%% Composition Type of species%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    compositiontype = 1; %%: partial densities [kg/m3] 

   % compositiontype = 2; %%: partial number densities [1/m3] 

   % compositiontype = 3; %%: molar fractions [-] 


%%%% Atmosphere Composition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%here is important to check the order of the species between the input
%solution file - the mixture file in mutation - the mixture in borat

    mixture_name = 'Mars_CO2_N_ionized';     %%: Gas mixture of case
    
    
%%%% State Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    state_model_str = 'ChemNonEq1T';    %%: 
    %state_model_str = 'ChemNonEqTTv';   %%: 


%%%% Link to Mutation (Linux) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% this is to find the mutation++ matlabinterface script called "retrieve_and_sort.m", where also mppcalc needs to be
    addpath('mutationlink');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Plot Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(groot,'DefaultFigureColormap',gray)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp_edensity=0;
disp_coll_freq=0;
disp_X=0;
disp_Z=0;
disp_mu=1;
disp_kappa=0;
disp_mesh=1;
writefig=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Display of most important input information  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%verification of input information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (dir1==dir2 && pooo1(1)==pooo2(1) && pooo1(2)==pooo2(2))
    error('dir1==dir2 and main points of origin are identical');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('-----------------------------------------------------------------------------------------');
fprintf('\nThis is the %s (%s) V%.2f',NAME,NAMESHORT,VERSION);
fprintf('\n-----------------------------------------------------------------------------------------');
fprintf('\n \t main variables name and values:');
fprintf('\n \t  flowfieldfilename=%s',flowfieldfilename);
fprintf('\n \t  scdir=%.0fdeg scdirr=%.0fdeg', scdir,scdir_range);
fprintf('\n \t  pooo1(1)=%.2fm pooo1(2)=%.2fm pooo2(1)=%.2fm pooo2(2)=%.2fm ',pooo1(1),pooo1(2),pooo2(1),pooo2(2));
fprintf('\n \t  dir1=%.0fdeg dir2=%.0fdeg',dir1,dir2);
fprintf('\n \t  stepsize=%.2f maxsteps=%d maxangles=%d',ss,maxsteps,maxangles);

pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Read flowfield; compute and add optical properties %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  analyse flowfield %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n %.2f s step 1 - reading flowfield variables, computing optical properties',cputime-initime)

if read == 1
    
    domain=readflowfield_tecplot(flowfieldfilename,f, collisionmodel,shrinkfactor,indexofelectronspecies,Species,pltformat,mixture_name,state_model_str,nb_species,indexof1stspecies,nb_temperatures,indexof1sttemperature,compositiontype,indexoftotaldensity);

elseif import ==1
    
    domain=import_flowfield_tecplot(flowfieldfilename,f, collisionmodel,shrinkfactor,indexofelectronspecies,Species,pltformat,mixture_name,state_model_str,nb_species,indexof1stspecies,nb_temperatures,indexof1sttemperature,compositiontype,indexoftotaldensity);

else
    
    error('flowfield reading procedure not identified')

end

domain.limitx1=limitx1;
domain.limitx2=limitx2;

domain.limity1=limity1;
domain.limity2=limity2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% compute important parameters

% compute max density in domain
max_edensity=0;

for z=1:domain.nozones %for each zone
    if max( domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:) ) > max_edensity
        max_edensity=max( domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:) );
    end
end

criticaldensity=f^2/80.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 1. Plot Domain Before Raytracing %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  plot flowfield %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if import ==1
    
    if flag_exportdomain==1
        
        export_flowfield_tecplot(domain,output_filename,collisionmodel,criticaldensity,max_edensity);
        command=strcat('mv *.szplt',32, runFolder); %32 is the space code
        unix(command);
    end
    
end

plotdomain(domain,criticaldensity,max_edensity,disp_edensity,disp_coll_freq,disp_mu,disp_kappa,disp_X,disp_Z,disp_mesh,writefig,runFolder);


pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display characteristics
fprintf('\n \t mesh characteristics:');
fprintf('\n \t  unstructured_mesh=%d  number_of_zones=%d',domain.unstructured,domain.nozones);
fprintf('\n \t  shrinkfactor %.1f',shrinkfactor);

for z=1:domain.nozones %for each zone
    fprintf('\n \t  zone=%d: N=%d E=%d surface_area=%.5fm^2',z,domain.( strcat('zone',num2str(z)) ).N,domain.( strcat('zone',num2str(z)) ).E,domain.( strcat('zone',num2str(z)) ).v);
end

fprintf('\n \t thermodynamic/optical characteristics:');
fprintf('\n \t  f %1.3e criticaldensity=%1.3e maxdensity=%1.3e',f,criticaldensity, max_edensity)
fprintf('\n \t  collisionmodel=%d',collisionmodel);
fprintf('\n %.2f s step 1 - reading flowfield variables - done',cputime-initime)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2. Ray-tracing: compute paths of rays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if solverSnell==1

    [itdir, itpo,symmetrylineencounter]=raytracing(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, absorptionlimits(2)-10,symmetryline,ss,f);

    fprintf('\n %.2f s step 2 - ray-tracing - Snell Law solver ',cputime-initime)

    
elseif solverEikonal==1

    [itdir, itpo,symmetrylineencounter]=eikonal2D(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, absorptionlimits(2)-10,symmetryline,ss,f);

    fprintf('\n %.2f s step 2 - ray-tracing - Eikonal solver ',cputime-initime)

else
    
    error('Rays integration solver not identified')
    
end
fprintf('\n %.2f s step 2 - ray-tracing - done',cputime-initime)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3. Draw ray-tracing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n %.2f s step 3 - plot ray-tracing',cputime-initime)

draw_raytracing(domain,plotscdir,scdir,scdir_range,collisionmodel,absorptionlimits,pooo1,pooo2,itpo,maxangles,itdir,symmetrylineencounter,writefig,runFolder)

fprintf('\n %.2f s step 3 - plot ray-tracing - done',cputime-initime)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%% 5. output final result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n %.2f s step 4 - output final results',cputime-initime)
%% 5.1 plot domain, spacecraft direction, optimum ray-tracing path, attenuation for that path
fprintf('\n');
%if symmetryline
%    final_dir_attenuation=zeros(2*maxangles,2); %% column one will contain final direction, column two will contain final attenuation
%else
final_dir_attenuation=zeros(maxangles,2); %% column one will contain final direction, column two will contain final attenuation
initial_dir_attenuation=zeros(maxangles,2);


for a=1:1:maxangles
    %% final attenuation for all rays
    if collisionmodel~=0 %% no attenuation if collisions/absorption was not accounted for
        final_dir_attenuation(a,2)=max(itpo(5, : , a ));
        final_dir_attenuation(a+1,2)=max(itpo(5, : , a ));
        
    else
        final_dir_attenuation(a,2)=0;
    end
    
    %end
    
    
    %% number of steps required to reach boundary, if reached
    %!add here more info
    %!  final ray density and initial angle density   @ initial angle
    %!	total antennuation per ray
    
    %% correct direction if reflected odd times at symmetry line %this for asymmetrical antenna, e.g. ExoMars
    if symmetryline
        if mod(symmetrylineencounter,2) %% odd number of reflections at symmetry line
            final_dir_attenuation(a,1)= -itdir(nnz(itdir(:,a)),a)+360;
        else
            final_dir_attenuation(a,1)= itdir(nnz(itdir(:,a)),a);
        end
    else
        final_dir_attenuation(a,1)= itdir(nnz(itdir(:,a)),a);
    end
    %% make sure that final direction is within range 0-360 degree
    if final_dir_attenuation(a,1)<0
        final_dir_attenuation(a,1)= final_dir_attenuation(a,1)+360;
    elseif final_dir_attenuation(a,1)>360
        final_dir_attenuation(a,1)=final_dir_attenuation(a,1)-360;
    else
        final_dir_attenuation(a,1)= final_dir_attenuation(a,1);
    end
    
    fprintf('\tn_steps: %d/%d \tinitial angle/finalangle %.1f/%.1f attenuation %.4f \n', nnz(itpo(1,:,a)),maxsteps,itdir(1,a), final_dir_attenuation(a,1), final_dir_attenuation(a,2) );
end

%% 5.2 plot attenuation over angles
%!this is not yet final, check dimensions and correctness
if collisionmodel~=0
    
    plot_attenuation(final_dir_attenuation,writefig,runFolder);
    
end
res = gcat(labindex);

cd (runFolder); %use the unix command? or it's fine like this
save(solutionFile);
cd ..;

command=strcat('mv *.dat',32, runFolder); %32 is the space code
unix(command);

%% end
fprintf('\n %.2f s step 5 - output final results - done',cputime-initime)

fprintf('\n Thank you for using %s. Exiting now.\n\n',NAMESHORT);
fprintf('-----------------------------------------------------------------------------------------\n');

%% --------------------------------------------------------------------------
save('raytracing_solution_KnappIRS_1mag.mat');%,'Y_int_1','Y_int_2','Y_int_3','YY','-append'); 
