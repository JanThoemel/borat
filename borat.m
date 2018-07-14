clear all; clc;format long;close all; initime=cputime; 
pause('off'); %there are pause points in the script, enable this option and they won't be active
%pause('on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NAME='BlackOut RAy-Tracer';NAMESHORT='BORAT';VERSION=1.2; 
% 23 Nov 2017
% by Jan Thoemel
%% system requirements:
%   -installed mutation++
%   -use on MS Windows if mutation dependent features (collision/absorption modeling) are not used
%   -tested on Matlab 2016b on linux 
%% release notes:
%   -
%% features:
%   -choice of 1) CoolFluid 2) Sahadeo (Tecplot written) plt format
%   -option to use symmetry line at y=0
%   -tested for CO2 and air (ExoMars and Stardust case)
%   -structured/unstructured meshes, tecplot format, automatic detection, some data format conventions apply
%   -use of Mutation++, tested for the air11 and CO2 (Stardust and ExoMars case)
%	-several initial equistiant points along a line between an initial and and final point
%	-equistant directions between and inital and a final direction 
%   -optical properties i.e. refractive index and absorption coefficients based on Davies
%       *w/o collisions-done
%       *collisions frequency from mutation
%   -triangle for interpolation is identified through insidedomain function
%   -multizones based on struct-arrays containing: 
%       a) flag structured/unstructured
%       b) no of zones
%       c) for each zone: number of points: i) i and j for structured ii) total no of points if unstructured
%       d) mesh points and values (x, y, and others in particular number densities of al species)
%       e) place holder for ri and cf
%   -basic ray-tracing analysis incl. attenuation and polar plot (not validated yet)
%   -critical density plotted
%% features to come
%   -count reflections on symmetry axis (for directivity) and total reflections, implement
%   symmetric, halfsymmetric(symmetric flow field, asymmetric antenna) polarplot
%   -check warnings
%   -use of alternative gradient method ( maybe: https://nl.mathworks.com/matlabcentral/fileexchange/36837-trigradient-m') for Snell's law
%   -get molar masses from mutation
%   -validation (,e.g. against Torino results)
%   -parametrization of snell's law gradient finding, i.e. angle (or use matlab built-in function)
%   -automatic detection of 1) CoolFluid 2) Sahadeo (Tecplot written) plt format
%	-reflection on vehicle surface
%   -further ray-tracing analysis/chapter 5
%   -dynamic stepsize for raytracer based on ri ratio po_0/po_far
%	-count number and output of total reflections because they change the phase
%% use comments:
%   -star dust air case works in 1 T mode (mutation) only
%   -compute regions that are above cut-off frequency and/or compute plasma frequency and/or critical electron density
%   -check for use of global variables, avoid this
%   -consider to store data in cell arrays or array struct instead of layered struct to improve performance. Currently the latter is done but
%   -optical properties implemented only for plasma and not for neutral gas/ should this be done?
%% short installation instructions
%   0) verify that you use linux and have a Matlab version 2016b or later using the Matlab command "ver"
%   1) place this Matlab script file in a convenient location
%   2) create the directory "cases"
%   3) install mutation++, see https://www.mutationpp.org/docs/html/installation.html
%   4) create the directory "mutationlink"
%   5) place "retrieve_and_sort.m" and "mppcalc.cpp" and the mppcalc "Makefile" into directory "mutationlink"
%   6) compile mppcalc.cpp using the makefile
%   7) place CFD files in directory "cases", make sure that format is sufficient
%   8) adapt all input parameters of BORT below to run the case, make sure they makes sense (obviously) and
%   in particular make sure that a) your settings b) the CFD input file c) the mutation settings are
%   compatible
%   9) now, you are ready to run this Matlab script
%   X) Remark: skip steps 4)-6) if you don't use mutation++, i.e. you are not accounting for absorption 
%   caused by collisions. Latter would be computed by mutation++. In this case you can run BORT also
%   on MS Windows
%% licence:
%    -no licence, use as you please, however:
%       1) no warranty whatsoever
%       2) refer to the source and the first paper (to be written)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% tgo direction in polar plot
%collision freq.label
% artifact close to vehicle surface


%%%%%% 0. definition of input data and information %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% basic verification and screen display thereof %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% transmission frequency, Hz
f=400.0e6;                  %% ExoMars UHF frequency (TBC)
%f=400.0e6; 				%% MSL UHF frequency
%f=8401.4e6;				%% MSL X-Band frequency, transmit of MGA (source: MSL Telecommunication System Design by Andre Makovsky et al.)
%f=7150.8e6;				%% MSL X-Band frequency, receive of MGA (source: MSL Telecommunication System Design by Andre Makovsky et al.)
%f=49.97e6;				%% BRAMS radar frequency
%f=1575.42e6; 				%% GPS L1 (wikipedia 20 may 2017)
%f=1227.60e6; 				%% GPS L1 (wikipedia 20 may 2017)
%f=1575.42e6; 				%% Galileo E1 (http://galileognss.eu/galileo-frequency-bands/ 20 may 2017)
%f=1191.795e6; 				%% Galileo E5 (http://galileognss.eu/galileo-frequency-bands/ 20 may 2017)
%f=1278.75e6; 				%% Galileo E6 (http://galileognss.eu/galileo-frequency-bands/ 20 may 2017)
%f=1618.85e6; 				%% Iridium lower end (wikipedia 20 may 2017)
%f=1626.5e6; 				%% Iridium upper end (wikipedia 20 may 2017)
%f=2267.0e6;                %% TDRS ( MODELLING OF ANTENNA RADIATION PATTERN OF A RE-ENTRY VEHICLE IN PRESENCE OF PLASMA G. Vecchi(’), M. Sabbadini(’1, R. Maggiora(’), A. Siciliano(’))

%% flowfield file related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flowfieldfilename='cases/th102blok_cat/th10.plt';                                    %%jan thesis development case
%flowfieldfilename='cases/stardust/Star.plt.LTE_FEF';                                 %%star dust

%flowfieldfilename='cases/ExoMars/t17s/ExoMars_t17s_LARS.dat';                         %% ExoMars case, format as sent by Sahadeo
% XX flowfieldfilename='cases/ExoMars/t17s2OCorN/ExoMars_t17s_LARS.dat';                         %% ExoMars case, format as sent by Sahadeo

%flowfieldfilename='cases/ExoMars/t22s2OCorN/ExoMars_t22s_LARS.dat';                         %% ExoMars case, format as sent by Sahadeo

%XXflowfieldfilename='cases/ExoMars/t27s2OCorN/ExoMars_t27_LARS.dat';                         %% ExoMars case, format as sent by Sahadeo

flowfieldfilename='cases/ExoMars/t38s2OCorN/ExoMars_t38_LARS.dat';                         %% ExoMars case, format as sent by Sahadeo

%XXflowfieldfilename='cases/ExoMars/t50s2OCorN/ExoMars_t50s_LARS.dat';                        %% ExoMars case, format as sent by Sahadeo

%flowfieldfilename='cases/ExoMars/t73s2OCorN/ExoMars_t73s_LARS_SecondOrder.dat';       %% ExoMars case, format as sent by Sahadeo, includes velocities, 1 temperature
%XXflowfieldfilename='cases/ExoMars/t73sTTV2O/ExoMars_t73s_TTv_LARS.dat';                        %% ExoMars case, format as sent by Sahadeo

%XX flowfieldfilename='cases/ExoMars/t80s2OCorN/ExoMars_80s_LARS.dat';                         %% ExoMars case, format as sent by Sahadeo

%XXflowfieldfilename='cases/ExoMars/t85s2OCorN/ExoMars_t85s_LARS.dat';                         %% ExoMars case, format as sent by Sahadeo

%XXflowfieldfilename='cases/ExoMars/t92s2OCorN/ExoMars_t92_LARS.dat';                         %% ExoMars case, format as sent by Sahadeo


%flowfieldfilename='cases/ExoMars/t104s/ExoMars_t104s_LARS.dat';                      %% ExoMars case, format as sent by Sahadeo
%flowfieldfilename='cases/ExoMars/t123s/ExoMars_t123s_allCO2.dat'                     %% ExoMars case, format as sent by Sahadeo
%flowfieldfilename='cases/meteor/MultiLarsenSol_CSV_mod.csv'                          %% Stefano's first meteor data sent, 10th Nov. 2017, format slightly modified

%% tecplot file format: 1: as provided by Coolfluid, 2: as provided by Sahadeo after interpolation of LARSEN results, see end of this MATLAB script for an example
%pltformat=2;
pltformat=2;
%% Is y=0 a symmetry line for the raytracing? This requires the domain to be above and adjacent to it.
%% Choose this option if you model an axi-symmetric case where the axis is at y=0
%symmetryline=1;
symmetryline=1;
%% shrinkfactor for the matlab boundary function to determine concave hull of each zone (do we need one per zone?)
%shrinkfactor=0.1;          %%jan thesis development case
shrinkfactor=0.999;         %%ExoMars&full star dust front shield

%% ray-tracing related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% rays are started from line with the start and end point as defined here into a direction as described further down
%% rays are started equi-distant between that point and equi-angled in between the min and max angles

%% rays are started from line with the start and end point
%% pooo1 shall be left and up to pooo2
%% if they all rays shall start from one point only then pooo1=pooo2
%pooo1=[0.3 0.001]; 		%%jan thesis development test cast
%pooo2=[0.3 0.001]; 		%%jan thesis development test cast
%pooo1=[0.1 0.1];           %%jan thesis development test cast
%pooo2=[0.5 0.1];           %%jan thesis development test cast
%pooo1=[0.2 0.001]; 		%%jan thesis development test cast
%pooo2=[0.4 0.001]; 		%%jan thesis development test cast
%pooo1=[0.308 1.0];			%%front shield exomars
%pooo2=[0.308 1.0];			%%front shield exomars
%pooo1=[0.365 0.325];		%%full stardust 
%pooo2=[0.365 0.325];		%%full stardust 
%pooo1=[0.5 1];             %%full stardust radar
%pooo2=[4.8 1];             %%full stardust radar
%pooo1=[0.89936 0.73545];    %%temporary exomars
%pooo2=[0.89936 0.73545];    %%temporary exomars
%pooo1=[1.269 0.3395];       %%2nd temporary exomars
%pooo2=[1.269 0.3395];       %%2nd temporary exomars
pooo1=[0.6987+0.01 0.9321+0.01];       %%3rd temporary exomars
pooo2=[0.6987+0.01 0.9321+0.01];       %%3rd temporary exomars
%pooo1=[1 7.5];       %% meteor
%pooo2=[2500 7.5];       %% meteor

%%  on angle convention:
%%  cone angle as in Morabito "The Mars Science Laboratory EDL
%%  Communications Brownout and Blackout at UHF" figure 9. i.e. from antiram over
%%  zenit (take care zenit is shown downwards in figure 9)

%% min and max ray directions, deg, dir1 shall be counterclockwise to dir2
%dir1=135;					%%development jan thesis case
%dir2=45;					%%development jan thesis case
%dir1=270;					%%development jan thesis case
%dir2=250;					%%development jan thesis case
%dir1=90;					%%development jan thesis case
%dir2=90;					%%development jan thesis case
%dir1=180;					%%exomars front shield
%dir2=135;					%%exomars front shield
dir1=43+80;					%%ExoMars full case TBC
dir2=43-80;					%%ExoMars full case TBC
%dir1=300;					%% meteor
%dir2=300;					%% meteor

%% plot space craft direction and range in refractive index+rays plot?
%plotscdir=0;
plotscdir=1;

%% spacecraft direction, cone angle, deg
%scdir=46.0443553950332; %t=-8
%scdir=42.7455125355700;	%t=17
%scdir=42.092770314134340;	%t=22
%scdir=41.4432946652926; %t=27
scdir=39.9682994082497; %t=38.50
%scdir=38.5434967075915; %t=50
%scdir=36.1398791092402; %t=73
%scdir=35.6479509374479; %t=80	
%scdir=35.38570638831026; %t=80	
%scdir=35.1559757351031; %t=92
%scdir=35.1304124788995; %t=104	
%scdir=35.8482606257494; %t=123	

%% spacecraft direction range(angular), deg
scdir_range=5;

%% max no of angles for raytracing, choose value that allows integer step size, i.e. 179,89,59 etc
%maxangles=15;
%maxangles=31;
%maxangles=15;
maxangles=25;
%% stepsize [m], the stepsize should be around 10 times greater then the CFD cell size (at relevant locations)
%ss=0.01;
ss=0.02;
%% max steps for ray-tracing, the product of stepsize (above) and maxsteps shall be large enough that the ray tracing reaches the boundary,
%% i.e.: ss*maxsteps>2*size of domain
maxsteps=2000;

%% physics modeling related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% collision model, for explanation see in-line in code in optical properties (0: no collisions, 1: mutation)
collisionmodel=1;

%% absorptionlimits for ploting; absorptionlimits(2)-10 is used to interrupt raytracing because it is assume the signal strength is too low to continue
absorptionlimits=[1e-6 40];

%% mutation settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% more mutation related settings are hard coded in function "readflowfield_tecplot"
%% this is to find the mutation++ matlabinterface script called "retrieve_and_sort.m", where also mppcalc needs to be
addpath('mutationlink');
%% make sure that order of species is identical in tecplot file and in mutation mixture file
%% description of inputfile
%nb_species=11;%meteor, stardust
nb_species=14;%% ExoMars
indexof1stspecies=3;
nb_temperatures=1;
%nb_temperatures=2;
%indexof1sttemperature=14; %meteor
%indexof1sttemperature=16;
indexof1sttemperature=19;% exomars 
%indexof1sttemperature=17;% exomars if velocities not provided

%not necessary to set if compositiontype = 2, see just below
indexoftotaldensity=20;

%% compositiontype: 1: partial densities [kg/m3] 2: partial number densities [1/m3] 
compositiontype=2;

%mixture_name = 'air11Jan';
mixture_name = 'Mars_CO2_N_ionized';
state_model_str = 'ChemNonEq1T';   % State model for chemical non-equilibrium single Temperature (For 2 temperatures use ChemNonEqTTv)
%state_model_str = 'ChemNonEqTTv';   % State model for chemical non-equilibrium single Temperature (For 2 temperatures use ChemNonEqTTv)

%% what plots shall be displayed & shall plotted figures be saved as eps files?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp_edensity=1; disp_coll_freq=1; disp_X=1; disp_Z=1; disp_mu=1; disp_kappa=1; disp_mesh=0;writefig=1;
%disp_edensity=0; disp_coll_freq=0; disp_X=0; disp_Z=0; disp_mu=0; disp_kappa=0; disp_mesh=0;writefig=0;
%disp_edensity=1; disp_coll_freq=0; disp_X=0; disp_Z=0; disp_mu=0; disp_kappa=0; disp_mesh=0;writefig=1;
%disp_edensity=1; disp_coll_freq=0; disp_X=0; disp_Z=0; disp_mu=0; disp_kappa=0; disp_mesh=0;writefig=1;

%% define colormap for surface plots
%set(groot,'DefaultFigureColormap',parula) 
set(groot,'DefaultFigureColormap',gray)



%% verification of input information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (dir1==dir2 && pooo1(1)==pooo2(1) && pooo1(2)==pooo2(2))
	error('dir1==dir2 and main points of origin are identical');
end


%% display of most important input information onto screen
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
%%%%% 1. read flowfield; compute and add optical properties %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  analyse flowfield, display flowfield  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n %.2f s step 1 - reading flowfield variables, computing optical properties',cputime-initime)

%% read domain
domain=readflowfield_tecplot(flowfieldfilename,f, collisionmodel,shrinkfactor,pltformat,mixture_name,state_model_str,nb_species,indexof1stspecies,nb_temperatures,indexof1sttemperature,compositiontype,indexoftotaldensity);

%% compute max density in domain
max_edensity=0;
for z=1:domain.nozones %for each zone
    if max(    domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)) > max_edensity
        max_edensity=max(    domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:));
    end
end
%% compute important parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
criticaldensity=f^2/80.5; 

%% plot domain before raytracing starts
plotdomain(domain,criticaldensity,max_edensity,disp_edensity,disp_coll_freq,disp_mu,disp_kappa,disp_X,disp_Z,disp_mesh,writefig);
pause

%% display characteristics
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

%%%%% 2.ray-tracing: compute paths of rays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n %.2f s step 2 - ray-tracing',cputime-initime)
[itdir, itpo,symmetrylineencounter]=raytracing(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, absorptionlimits(2)-10,symmetryline,ss);
fprintf('\n %.2f s step 2 - ray-tracing - done',cputime-initime)


%%%%% 3.draw ray-tracing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n %.2f s step 3 - plot ray-tracing',cputime-initime) 
%% compute direction and range to space craft
if plotscdir
    if ( pooo1(1)==pooo2(1) && pooo1(2)==pooo2(2) ) %% compute and show spacecraft direction and range only if raytracing start from antenna location i.e. pooo1=pooo2
        for it=1:3
          scdir_x(it)=pooo1(1)+cos(scdir/180*pi)*(it-1)*3;                 scdir_y(it)=pooo1(2)+sin(scdir/180*pi)*(it-1)*3;
          scdirr1_x(it)=pooo1(1)+cos((scdir+scdir_range)/180*pi)*(it-1)*3; scdirr1_y(it)=pooo1(2)+sin((scdir+scdir_range)/180*pi)*(it-1)*3;
          scdirr2_x(it)=pooo1(1)+cos((scdir-scdir_range)/180*pi)*(it-1)*3;	scdirr2_y(it)=pooo1(2)+sin((scdir-scdir_range)/180*pi)*(it-1)*3;
        end
    end
end
%% plot rays and spacecraft direction/direction range within refractive index contour plot
figure       
		hold on;
    	V=[1e-17,1e-16,1e-15,1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-3:0.02:1-1e-3,1-1e-3,1-1e-4,1-1e-5,1-1e-6,1-1e-7,1-1e-8,1-1e-9,1-1e-10,1-1e-11,1-1e-12,1-1e-13,1-1e-14,1-1e-15,1-1e-16,1-1e-17,1];
        %V=[0:0.1:0.9,0.9:0.01:1];
		%V=[0.09:0.01:0.99];
        %axis equal;
        for z=1:domain.nozones
			%%plot flowfield with ri, raytraces and s/c direction/range
            
            %%plot ri contour plot
            if 0 
				%% plot boundary
				plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)                           
				%% plot ri
                [C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),V);
            end  
            %%plot ri surface plot
			if 1 
				%% plot boundary
				plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)                           
				%% plot ri
				h=trisurf(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)));            
				set(h, 'edgecolor','none');
			end
        end
         
        z_max = max(max(get(h,'Zdata')));
        scdir_z=[z_max,z_max,z_max];       

        %%plot rays in domain
        itpoz=z_max*ones(   size( itpo(1,:,1),2)   ,1);%z coordinates of lines to be plotted
        for a=1:maxangles                     
            %!can they be colored according to attenuation?
            plot3(  itpo(1,1:nnz(itpo(1,:,a)),a)  ,  itpo(2,1:nnz(itpo(1,:,a)),a),itpoz(1:nnz(itpo(1,:,a))),'k-','LineWidth',1,'color',[105/256 105/256 105/256] );
        end
        %%plot s/c direction and direction range
        if plotscdir
            if ( pooo1(1)==pooo2(1) && pooo1(2)==pooo2(2) )%% compute and show spacecraft direction and range only if raytracing start from antenna location i.e. pooo1=pooo2
                plot3(scdir_x,scdir_y,scdir_z,'k-.','color',[20/256 20/256 20/256]);
                %plot3(scdirr1_x,scdirr1_y,scdir_z,'k--','color',[20/256 20/256 20/256]);
                %plot3(scdirr2_x,scdirr2_y,scdir_z,'k--','color',[20/256 20/256 20/256]);
            end
        end
        text(1.5,3.1,z_max+1,'direction to TGO');
        set(h, 'edgecolor','none');
        xlim([-0.17 6]);ylim([0 5.185]);caxis([0 1])

        hcb=colorbar;title(hcb,'\mu');
        box off;grid off;
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])       
        hold off;
if writefig
    savefig('raytracing')
    print('-painters','raytracing','-depsc')            
    print('-painters','raytracing','-dpng')            
end

fprintf('\n %.2f s step 3 - plot ray-tracing - done',cputime-initime) 


%%%%% 4. ray-tracing analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n %.2f s step 4 - ray-tracing analysis',cputime-initime) 
%% plot for ri over path length
figure       
    hold on
    for a=1:maxangles
        plot(itpo(3,1:nnz(itpo(1,:,a)),a),itpo(4,1:nnz(itpo(1,:,a)),a))
    end
    xlabel('path length [m]');ylabel('refractive index [-]');axis([0 inf -0.05 1.05]);
    hold off
    if writefig
        savefig('path_mu')        
        print('-painters','path_mu','-depsc')            
        print('-painters','path_mu','-dpng')            
    end
%% plot attenuation over path length    
if collisionmodel~=0 %% no plot if collisions/absorption was not accounted for
    figure
        for a=1:maxangles
            semilogy(  itpo(3,1:nnz(itpo(1,:,a)),a)  ,  itpo(5,1:nnz(itpo(1,:,a)),a)  );
            hold on
        end
        xlabel('path length [m]');ylabel('attenutation [dB]');axis([0 inf absorptionlimits]);
        %% absorption cut off
        semilogy([0 inf], [absorptionlimits(2)-10 absorptionlimits(2)-10]);
        if writefig
            savefig('path_kappa')            
            print('-painters','path_kappa','-depsc')            
            print('-painters','path_kappa','-dpng')            
        end
        hold off
end

fprintf('\n %.2f s step 4 - ray-tracing analysis - done',cputime-initime)


%%%%% 5. output final result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n %.2f s step 5 - output final results',cputime-initime) 
%% 5.1 plot domain, spacecraft direction, optimum ray-tracing path, attenuation for that path
fprintf('\n');
%if symmetryline
%    final_dir_attenuation=zeros(2*maxangles,2); %% column one will contain final direction, column two will contain final attenuation
%else
    final_dir_attenuation=zeros(maxangles,2); %% column one will contain final direction, column two will contain final attenuation
%end
%if symmetryline
%    for a=1:2:2*maxangles
%        %% final attenuation for all rays
%        if collisionmodel~=0 %% no attenuation if collisions/absorption was not accounted for
%            final_dir_attenuation(a,2)=max(itpo(5, : , a ));
%            final_dir_attenuation(a+1,2)=max(itpo(5, : , a ));
%        else
 %           final_dir_attenuation(a,2)=0;
%        end
%    end
%else
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
    
    fprintf('\tnos: %d/%d \tinitial angle/finalangle %.1f/%.1f attenuation %.4f \n', nnz(itpo(1,:,a)),maxsteps,itdir(1,a), final_dir_attenuation(a,1), final_dir_attenuation(a,2) );
end
%% 5.2 plot attenuation over angles
%!this is not yet final, check dimensions and correctness
if collisionmodel~=0
    final_dir_attenuation=sortrows(final_dir_attenuation);
    figure
        polarplot( [final_dir_attenuation(:,1)'*pi/180 final_dir_attenuation(1,1)*pi/180], [final_dir_attenuation(:,2)' final_dir_attenuation(1,2)] );
        hold on
        rlim([1e-5 1e-1])
        %polar([scdir scdir] [rlim(1) rlim(2)],'k-.','color',[20/256 20/256 20/256]);
        hold off
        if writefig
            savefig('polar_mu')
            print('-painters','polar_mu','-depsc')            
            print('-painters','polar_mu','-dpng')            
        end
        
end
fprintf('\n %.2f s step 5 - output final results - done',cputime-initime)

fprintf('\n Thank you for using %s. Exiting now.\n\n',NAMESHORT);
fprintf('-----------------------------------------------------------------------------------------\n');

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------------------- LEVEL 1 ---------------------------------
%------------------------------- FUNCTIONS---------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function domain=readflowfield_tecplot(flowfieldfilename,f,collisionmodel, shrinkfactor,pltformat,mixture_name,state_model_str,nb_species,indexof1stspecies,nb_temperatures,indexof1sttemperature,compositiontype,indexoftotaldensity)
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
    options.erase_temporary_file = true;
    %options.erase_temporary_file = false;
    %options.display = true;
    options.display = false;
    partialdensityfix=1e-30;    %if any partial density is below this value, set it to this value
    temperaturefix=150;         %if any temperature is below this value, set it to this value
    
    %% molar mass of each species
    MM=zeros(nb_species,1);   
    if  strcmp(mixture_name,'air11Jan')
    %% air11
    	MM(1)=5.48579909070e-7; MM(2)=7e-3; MM(3)=8e-3; MM(4)=14e-3; MM(5)=15e-3; MM(6)=16e-3; MM(7)=14e-3;	MM(8)=15e-3;	MM(9)=14e-3;	MM(10)=16e-3;	MM(11)=8E-3;
        %em                      N           O           N2           NO           O2           N2p              NOp             Np              O2p             Op
    elseif   strcmp(mixture_name ,'Mars_CO2_N_ionized')
    %% Mars_CO2_N
        MM(1)=5.48579909070e-7; MM(2)=44e-3; MM(3)=28e-3; MM(4)=12e-3; MM(5)=14e-3; MM(6)=16e-3; MM(7)=32e-3;	MM(8)=28e-3;	MM(9)=30e-3;	MM(10)=28e-3;	MM(11)=30E-3;    MM(12)=12e-3;	MM(13)=16e-3;	MM(14)=32E-3;
        %e-                      CO2         N2           C            N            O            O2             CO              NO              CO+             NO+             C+              O+              O2+
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
                        domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i)*domain.( strcat('zone',num2str(z)) ).variables(indexoftotaldensity,i)*nA/molarmasselectron;
                   elseif compositiontype==2
                        domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,i)  =  domain.( strcat('zone',num2str(z)) ).variables(indexof1stspecies,i);
                   else
                       error('unknown compositiontype: %d',compositiontype);
                   end
                   %%collisionfrequency placeholder
                   domain.( strcat('zone',num2str(z)) ).variables(domain.nova+2,i)  = 0;
            end

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

            for i=1:domain.( strcat('zone',num2str(z)) ).N
                   %% optical properties, i.e. refractive index and absorption coefficient
                   [domain.(strcat('zone',num2str(z))).variables(domain.nova+3,i),domain.(strcat('zone',num2str(z))).variables(domain.nova+4,i),domain.(strcat('zone',num2str(z))).variables(domain.nova+5,i),domain.(strcat('zone',num2str(z))).variables(domain.nova+6,i)]  =  opticalproperties(f,collisionmodel, domain.(strcat('zone',num2str(z))).variables(:,i) ); 
            end
            pause
            %% read elements (this necessary to proceed to the next zone)
            domain.( strcat('zone',num2str(z)) ).elements=zeros(domain.( strcat('zone',num2str(z)) ).E,4 );
            for i=1:domain.( strcat('zone',num2str(z)) ).E
                   line=fgetl(id);
                   domain.( strcat('zone',num2str(z)) ).elements(i,:)=sscanf(line,'  %d');
            end
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
            dimline=fgetl(id);
        end % while eof
    fclose(id);
    domain.nozones=z;
    domain.nova=domain.nova+6; %%make this automatic ???
end %% function read_flowfield
 

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [itdir, itpo,symmetrylineencounter]=raytracing(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, absorptioncutoff,symmetryline,ss)
%% this function traces rays 

%%	input
%%
%%	output
%%
%%  internal parameters
%%  iteration message
    msg=' ';
    symmetrylineencounter=zeros(maxangles,1);
    
    %% trace from around wall normal
    if maxangles==1
        anglestep=0;
    else
        anglestep=(dir1-dir2)/(maxangles-1);
		x_step=(pooo1(1)-pooo2(1))/(maxangles-1);
		y_step=(pooo1(2)-pooo2(2))/(maxangles-1);
    end
    %%as alternative, angles for ray-tracing from wall to wall, deg
    %anglestep=180/(maxangles+1);
    %% iteration points: 1: x, 2: y, 3: pathlength, 4: ri, 5: ec
    itpo=zeros(5,maxsteps,maxangles);
    %% iteration directions
    itdir=zeros(maxsteps,maxangles);
    %%starting directions and points for all rays
    if maxangles==1
        itdir(1,1)=dir1;
		itpo(1,1,1)=pooo1(1);
		itpo(2,1,1)=pooo1(2);
        %% starting point refractive index for all rays
        itpo(4,1,:)=interpolation(domain,itpo(:,1,1),domain.nova-1);
        %% starting point absorption coefficient for all rays
        itpo(5,1,:)=interpolation(domain,itpo(:,1,1),domain.nova);
    else
        for a=1:maxangles
            %% starting direction (fan) 
            itdir(1,a)=dir1-anglestep*(a-1);
        	%% initialpoints
			itpo(1,1,a)=pooo1(1)-x_step*(a-1);
			itpo(2,1,a)=pooo1(2)-y_step*(a-1);
            %% starting point refractive index for all rays
            itpo(4,1,a)=interpolation(domain,itpo(:,1,a),domain.nova-1);
            %% starting point absorption coefficient for all rays
            itpo(5,1,a)=interpolation(domain,itpo(:,1,a),domain.nova);
        end %% for
    end %%if

    for a=1:maxangles
        s=1;msg=' ';
        %% while loop to check if ray-tracing stepping leaves 1) domain 2) maxsteps are exceeded or attenuation is above threshold
        %% if so go to next ray
		%! evaluate and improve efficiency of this loop

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        while checkifinsidedomain(itpo(1:2,s,a)',domain)==1 && s<maxsteps  && itpo(5,s,a)<absorptioncutoff
                    s=s+1;     
                    po_far(1)=itpo(1,s-1,a)+cos(itdir(s-1,a)/180*pi)*1.5*ss;           
                    po_far(2)=itpo(2,s-1,a)+sin(itdir(s-1,a)/180*pi)*1.5*ss; 
%                    if not(checkifinsidedomain(po_far,domain)) 
                    %! check this logic!
                    if not(  checkifinsidedomain(po_far,domain)  ) && not(po_far(2) < 0 && symmetryline)
                         fprintf('\tpo_far outside')                        
                         break;
                    end
                    %if po_far outside any zone - break raytracing of this ray     

                    po_temp(1)=itpo(1,s-1,a)+cos(itdir(s-1,a)/180*pi)*ss;           
                    po_temp(2)=itpo(2,s-1,a)+sin(itdir(s-1,a)/180*pi)*ss;
%                    if not(checkifinsidedomain(po_temp,domain))
                    %! check this logic!
                    if not(  checkifinsidedomain(po_temp,domain)  ) && not(po_temp(2) < 0 && symmetryline)                        
                        fprintf('\t po_temp outside');
                        break;
                    end

                    %% if symmetryline is on 
                    %% check if this the right way to implement the symmetry line
                       if po_temp(2)<0
                           itpo(1,s,a)=itpo(1,s-1,a);
                           itpo(2,s,a)=itpo(2,s-1,a);
                           itpo(3,s,a)=itpo(3,s-1,a)+sqrt(  (itpo(1,s,a)-itpo(1,s-1,a))^2 + (itpo(2,s,a)-itpo(2,s-1,a))^2   );
                           itpo(4,s,a)=interpolation(domain,itpo(1:2,s,a)',domain.nova-1);
                           %% integrate absorption coefficients
                           itpo(5,s,a)= itpo(5,s-1,a)+interpolation(domain,itpo(1:2,s,a)',domain.nova)*(itpo(3,s,a)-itpo(3,s-1,a)); 
                           itdir(s,a)=-itdir(s-1,a);
                           symmetrylineencounter(a)=symmetrylineencounter(a)+1;
                           continue
                       end

                    %!if po_far in next zone and po_temp in current zone, shorten stepsize
                    %!if po_far and po_temp in next zone enlarge step size

                    ri_po_temp=interpolation(domain,po_temp,domain.nova-1);

                    po_snellbase(1)=itpo(1,s-1,a)+cos(itdir(s-1,a)/180*pi)*ss/2;           
                    po_snellbase(2)=itpo(2,s-1,a)+sin(itdir(s-1,a)/180*pi)*ss/2; 

                    for it_del=1:size(msg,2)
                        fprintf('\b');
                    end
                    msg=sprintf('\n \ts=%d/%d a=%d/%d inidir=%.0f/(%.0f - %0.f)',s,maxsteps,a,maxangles,itdir(1,a),itdir(1,1),itdir(1,maxangles));
                    fprintf(msg);
 
                    %% call snellslaw
                    itdir(s,a)=snellslaw(po_snellbase,itdir(s-1,a),ss, domain,itpo(4,s-1,a), ri_po_temp);

                    itpo(1,s,a)=po_snellbase(1)+cos(itdir(s,a)/180*pi)*ss/2;           
                    itpo(2,s,a)=po_snellbase(2)+sin(itdir(s,a)/180*pi)*ss/2; 
                    %{          
                    %% use the following for improved robustness
                    if not(checkifinsidedomain(itpo(1:2,s,a)',variables))
                        fprintf('\titpo outside')
                        break;
                    end
                    %}
                    %% compute length of path
                    itpo(3,s,a)=itpo(3,s-1,a)+sqrt(  (itpo(1,s,a)-itpo(1,s-1,a))^2 + (itpo(2,s,a)-itpo(2,s-1,a))^2   );
                    %% save ri over path
                    itpo(4,s,a)=interpolation(domain,itpo(1:2,s,a)',domain.nova-1);
                    %% integrate absorption coefficients over path
                    itpo(5,s,a)=itpo(5,s-1,a)+interpolation(domain,itpo(1:2,s,a)',domain.nova)*(itpo(3,s,a)-itpo(3,s-1,a));
                    
        end % while loop over steps       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        if s==maxsteps          
            %% ray-tracing has not reached boundary
            fprintf('\t ray-tracing -  maximum steps reached');
        end % if
        if itpo(5,s,a)>=absorptioncutoff          
            %% while loop ended because absorptioncutoff was reached
            fprintf('\t attenuation %.3f above absorptioncutoff %.3f ',itpo(5,s,a),absorptioncutoff);
        end %if
        fprintf('\n');        
    end %for loop over angles
    %pause
end %% function raytracting

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function plotdomain(domain,criticaldensity,max_edensity,disp_edensity,disp_coll_freq,disp_mu,disp_kappa,disp_X,disp_Z,disp_mesh,writefig)
%% this function plot values in/of domain

%%  input
%%  domain
%%  disp_edensity,disp_coll_freq,disp_mu,disp_kappa,disp_X,disp_Z,disp_mesh
%%  output
%%  none
%%  internal parameters    
    nb_lines=10;
    
    if   disp_edensity  %% plot electrondensity/critical density                              
        figure
        hold on
       
        VD=[criticaldensity criticaldensity];
        
        
        for z=1:domain.nozones       
                %% plot boundary
                plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)            
                %% plot electrondensity countours
%                [C,h]=tricontour(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  ),nb_lines);                               
                %[C,h]=tricontour(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  ),VE);                               
                %[C,h]=tricontour(    domain.( strcat('zone',num2str(1)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(1)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(1)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(1)) ).variables(domain.nova-5,:)  ),VE);                               
                %VE=[1e13,1e14,1e15,1e16,1e17,1e18];
                
                h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  ) );            
                xlabel('X');ylabel('Y');zlabel('Z');
                set(h, 'edgecolor','none');
%                xlim([-0.17 6]);ylim([0 5.185]);caxis([1e16 1e17]);
                xlim([-0.35 6.1]);ylim([0 5.2]);caxis([1e16 1e17]);
                
                %h.Fill='on'
                %% plot critical electrondensity, Davis page 61, eq. 2.54 and k of page 71                  
                %{
                if max_edensity > criticaldensity
                    [D,i]=tricontour(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  ));
                    i(1).LineWidth=2;
                    i(1).EdgeColor='k';
                    %text(0,5,strcat('critical density=',num2str(criticaldensity,'%1.1e')));
                end
                %}
                hcb=colorbar;title(hcb,'X_e [1/m^3]          ');
                box off;grid off;
                set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])       

        end
        hold off
        if writefig
            savefig('e_density')            
            print('-painters','e_density','-depsc')            
            print('-painters','e_density','-dpng')            
        end
    end

    if disp_coll_freq    %%plot collision surface plot
        figure
        hold on

        for z=1:domain.nozones       
                %%plot boundary
                plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)                   
                %% plot coll freq
                %h=trisurf(domain.(strcat('zone',num2str(z))).delaunay,squeeze(domain.(strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)));            
                %  [CS,h]=tricontf(squeeze(domain.(strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,domain.(strcat('zone',num2str(z))).delaunay,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)));            
                %[CS,h]=TRICONTF(X,Y,M,Z)
                %h=trisurf(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)));            
                %[C,h]=tricontour(domain.(strcat('zone',num2str(z))).delaunay,squeeze(domain.(strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)),nb_lines);            
                %VC=[0,0.0001,0.001,0.01,0.10,1];

                h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)  ) );            
                set(h, 'edgecolor','none');
                xlim([-0.35 6.1]);ylim([0 5.2]);caxis([8000 12000]);

                hcb=colorbar;title(hcb,'collision freq. [1/s]                 ');
                box off;grid off;
                set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])       

        end
        hold off
        if writefig
            savefig('coll_freq')            
            print('-painters','coll_freq','-depsc')            
            print('-painters','coll_freq','-dpng')            
        end
    end 
    
    if disp_X    %% plot X
        figure
        hold on
        VX0=[1 1];
        for z=1:domain.nozones       
                %% plot boundary
                plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)            
                %% plot X
                %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)),nb_lines);
                %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)),VX);
                %VX=[1e-17,1e-15,1e-13,1e-11,1e-10,1e-8,1e-6,1e-4,1e-2,1e-1,2e-1:0.02:1-2e-1,1-1e-1,1-1e-3,1,3,5,7,9,20,40,60,80,100];
                
                h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)  ) );            
                set(h, 'edgecolor','none');
                xlim([-0.35 6.1]);ylim([0 5.2]);caxis([0 10])
                
                hcb=colorbar;title(hcb,'X [-]');
                box off;grid off;
                set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])       
                if max_edensity > criticaldensity
                    [D,i]=tricontour(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)  ),VX0);
                    i(1).LineWidth=2;
                    i(1).EdgeColor='k';
                    %text(0,4.5,'X=1');
                end

        end
        hold off
        if writefig
            savefig('X')            
            print('-painters','X','-depsc')            
            print('-painters','X','-dpng')            
        end
    end
    
    if disp_Z    %% plot Z
        figure
        %axis equal;
        hold on
        
        for z=1:domain.nozones       
                %% plot boundary
                plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)            
                %% plot Z
                %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-2,:)),nb_lines);
                %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-2,:)),VZ);
                %VZ=[1e-17,1e-16,1e-15,1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6];
                
                h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-2,:)  ) );                            
                set(h, 'edgecolor','none');
                xlim([-0.35 6.1]);ylim([0 5.2]);caxis([0 1e-6])
                
                hcb=colorbar;title(hcb,'Z [-]       ');
                box off;grid off;
                set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])       
        end
        hold off
        if writefig
            savefig('Z')            
            print('-painters','Z','-depsc')            
            print('-painters','Z','-dpng')            
        end
        
    end 
    
    
    if disp_mu    %% plot refractive index
        figure
        hold on
        for z=1:domain.nozones       
                %% plot boundary
                plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)            
                %% plot ri
                %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),nb_lines);
                %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),VRI);
                VRI=[1e-17,1e-15,1e-13,1e-11,1e-9,1e-7,1e-5,1e-3,1e-1,2e-1:0.02:1-2e-2,1-1e-1,1-1e-3,1-1e-5,1-1e-7,1-1e-9,1-1e-11,1-1e-13,1-1e-15,1-1e-17,1];
                
                h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)  ) );                            
                set(h, 'edgecolor','none');
                xlim([-0.35 6.1]);ylim([0 5.2]);caxis([0 1])
                
                hcb=colorbar;title(hcb,'\mu [-]');
                box off;grid off;
                set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])       
        end
        %img=imread('cases/ExoMars/ExoMarsHalfCapsule.png');
        %imagesc([1 1],[5 5],img);        
        hold off
        if writefig
            savefig('mu')            
            print('-painters','mu','-depsc')            
            print('-painters','mu','-dpng')            
        end
    end
    
    if disp_kappa    %% plot absorption coefficient
     
        figure
        hold on
        for z=1:domain.nozones       
                %% plot boundary
                plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)            
                %% plot kappa
                %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),nb_lines);
                %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),VRI);
                VRI=[1e-17,1e-15,1e-13,1e-11,1e-9,1e-7,1e-5,1e-3,1e-1,2e-1:0.02:1-2e-2,1-1e-1,1-1e-3,1-1e-5,1-1e-7,1-1e-9,1-1e-11,1-1e-13,1-1e-15,1-1e-17,1];
                
                h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova,:)  ) );                            
                set(h, 'edgecolor','none');
                xlim([-0.35 6.1]);ylim([0 5.2]);caxis([0 100]);
                
                hcb=colorbar;title(hcb,'\kappa [1/m]');
                box off;grid off;
                set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])       
        end
        %img=imread('cases/ExoMars/ExoMarsHalfCapsule.png');
        %imagesc([1 1],[5 5],img);        
        hold off
        if writefig
            savefig('kappa')
            print('-painters','kappa','-depsc')            
            print('-painters','kappa','-dpng')            
        end
    end
    
    
    if disp_mesh    %%plot mesh
        figure
        axis equal;
        hold on
        for z=1:domain.nozones       
                trimesh( domain.( strcat('zone',num2str(z)) ).delaunay  ,  domain.( strcat('zone',num2str(z)) ).variables(1,:)'   ,   domain.( strcat('zone',num2str(z)) ).variables(2,:)'   );
        end
        hold off
        if writefig
            savefig('mesh')            
            print('-painters','mesh','-depsc')            
            print('-painters','mesh','-dpng')            
        end        
    end
        
end % function plotdomain

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%------------------------------ LEVEL 2/3/4 -------------------------------
%------------------------------- FUNCTIONS---------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [X,Z,ri,ac]=opticalproperties(f,collisionmodel,variables)
%%  this function computes the:
%%  1) refractive index
%%  2) absorption coefficient, different models can be chosen or collision are neglected (then absorbtion coefficient is usually 0 except 1-X<0)

%%  input
%%  double  f               radiofrequency, Hz
%%  int     collisionmodel  parameter to chose whether collisions are taken into account and how they are  modeled.  
%%  struct  domain          all variables of the domain
%%  output
%%  struct  domain          like the input but with two additional variables (ri,ec)
%%  internal parameters
	%% Boltzmann constant, Google, m2 kg s-2 K-1
	kB=1.38064852e-23;
	%% Avogadro constant, wikipedia, mol−1
	nA=6.022140857e23;
	%% speed of light, wikipedia 8/7/2017, m/s
	c=299792458;
	
	%% 0. without accounting for collisions
    %% equation (2.78) (bookpage 71/pdf-page 89) of Ionospheric Radio Propagation by Kenneth Davies
    %% kk=e^2/(4 pi^2 e0 m)  
    kk=80.5;
    %% 1. with accounting for collisions
    %% equation (2.93) resolved for mu and k (bookpage 80/pdf-page 98) of Ionospheric Radio Propagation by Kenneth Davis
    %% this equation is to be preferred. It avoids negative mu 
    %% for definition of X,Z see Davies

   
    %! Do we need this switch? case 0 and 1 collapse if Z=0
    switch collisionmodel
        case 0 %% no collisions
            %% compute electrondensity
            electronnumberdensity=variables(size(variables,1)-5);
            X=kk*electronnumberdensity/f^2;
            %%compute optical properties
            if 1-X>=0
                %%refractive index
                ri=sqrt(1-X);
                %%absorption coefficient
                ac=0;
            elseif 1-X<0
                warning('ri_ac_domain_100')
                %%refractive index
                ri=0;
                %%absorption coefficient, check sign convention
                ac=2*pi*f/c*sqrt(X-1);                          
            end
        case 1 %% mutation collisions
            %% compute electrondensity
            electronnumberdensity= variables(size(variables,1)-5);
            %% compute collisionfrequency
            collisionfrequency=variables(size(variables,1)-4);
            %% compute X
            X=kk*electronnumberdensity/f^2;
            %% compute Z
            Z=collisionfrequency/(2*pi*f);
%{
            %% refractive index
            ri=            sqrt(     0.5*(    1-X/(1+Z^2)  +  sqrt(1 - 2*X/(1+Z^2) + 1/(1+Z^2)*X^2)   )       );
            %% absorption coefficient, neper/meter, 1 neper=8.69 dB, see Keneth Davies page 81
            ac=2*pi*f/c*  sqrt(      0.5*(   -1+X/(1+Z^2)  +  sqrt(1 - 2*X/(1+Z^2) + 1/(1+Z^2)*X^2)   )       );
%}
  
            root1=sqrt(1 - 2*X/(1+Z^2) + X^2/(1+Z^2));
            summand1=1-X/(1+Z^2);
            
            if abs(summand1)>abs(root1)
                %fprintf('\n sumcoorrection: X %1.30e Z %1.30e',X,Z);
                %fprintf('\n sumcoorrection: summand1 %1.30e rt1 %1.30e\n',summand1,rt1);
                %! why is the summand sometimes greater than the root, necessitating this fix?
                summand1=summand1*0.999999;
                %pause
            end
            ri=         sqrt(0.5* (summand1+root1));
            ac=2*pi*f/c*sqrt(0.5* (-summand1+root1));
            if  ~isreal(ri) || ~isreal(ac) 
                fprintf('\nX %f Z %f ri %f ac %f root1 %f summand1 %f',X,Z,ri, ac,root1,summand1);
                pause
            end
            %% absorption coefficient in dB, db
            ac=ac*8.69;
            %fprintf('\n X %e Z %e rt %3.10e sud %3.10e dif %3.10e ri %e ac %e',X,Z,rt1, summand1,(-summand1+rt1),ri,ac);
            %fprintf('\n X %e Z %e ',X,Z);
            %fprintf('\n  rt %3.30e sud %3.30e dif %3.30e ',rt1, summand1,(-summand1+rt1));
             %fprintf('\n  ri %e ac %e',ri,ac);
            %if not(isreal(ac))
            %   fprintf('\n ');
            %    ri
            %    ac
            %    rt1
            %    summand1
            %    pause
            %end
        case 2 %%  fake, i.e. fake electron density and fake or constant collisions
            %0. auxiliary variable for fake data
            %http://tutorial.math.lamar.edu/Classes/Alg/Parabolas.aspx
             %% auxiliary variables to generate fake data
            %maxvalue=max(max(domain.( strcat('zone',num2str(1)) ).variables(3,:,:)));
            h=0.15;k=1.7e15;a=-3.5e16;offset_x=0.35;offset_y=-0.01;
            radius=sqrt(   (variables(1)-offset_x)^2   +   (variables(2)-offset_y)^2    );
            %% compute electrondensity
            electronnumberdensity = a*(radius-h)^2+k  ;   
            %%compute collisionfrequency
            collisionfrequency=40000;
            %collisionfrequency=sqrt(a*(radius-h)^2+k );
            %compute X
            X=kk*electronnumberdensity/f^2;
            %compute Z
            Z=collisionfrequency/(2*pi*f);
            %%refractive index
            ri=0.5*(   1-X/(1+Z^2)  +  sqrt(  (1-X/(1+Z^2))^2 +Z^2*X^2/(1+Z^2)^2   ));
            %%absorption coefficient, neper/meter, 1 neper=8.69 dB, see Keneth Davies page 81, equation 2.98
            ac=2*pi*f/c/ri/2*Z*X/(1+Z^2);
            %%absorption coefficient, db
            ac=ac*8.69;
%        case 3 %% according to the Takahashi "Analyis of Radio Frequency Blackout for a Blunt-Body Capsule in Atmospheric Reentry Missions" (page 8, equation 14)
%            %%compute collisionfrequency
%            collisionfrequency=0;
%            for sp=1:nb_species-1 %sum over species but not over electrons
%                crosssection(sp)=dcs(sp,1)+dcs(sp,2);
%                collisionfrequency=collisionfrequency+numberdensity(sp)*pi*crosssection(sp)*sqrt(8*kB*T/pi/MM(1)/nA);
%            end
%            %compute Z
%            Z=collisionfrequency/(2*pi*f);
        otherwise
            error('unknown collisionmodel: %d',collisionmodel)
    end %% switch collision frequency model
               
end %% function optical properties


function snell_d2=snellslaw(snell_po_base,snell_d1,ss,domain,snell_ri1,snell_ri2)
%% 	application of Snell's law
%%	snellslaw: sin a1/sin a2=ri2/ri1
%% 	comment: no negative refraction considered, which is unlikely necessary

%% input
%%
%% output
%%
%% parameters
    totalreflection=0;

    %% find angle of isosurface
    surfaceangle=isolinefinding(domain,snell_po_base,ss);

    %% determine incoming angle, initial outgoing angle variable
    %! does this work for all situation e.g. downward rays? 
    if snell_d1>360
       %fprintf('\t snell_d1 > 360 deg, snell_d1: %.2f\n',snell_d1)
	   snell_d2=snell_d1-360;
    end
    if snell_d1<0
       %fprintf('\tsnell_d1 < 0 deg, snell_d1: %.2f\n',snell_d1)
	   snell_d1=snell_d1+360;
    end
    down=0;
    if snell_d1>180
        down=1;
        snell_d1=snell_d1-180;
    end
    
    if surfaceangle<snell_d1
        snell_a1=surfaceangle+270-(180+snell_d1);
    else
        snell_a1=surfaceangle+90-(180+snell_d1);
    end     
    
    if snell_a1>90 || snell_a1<-90
		 fprintf('\t snell_a2>90 deg or snell_a<-90, snell_a2: %.2f \n',snell_a1)
    end %% if
            
    snell_a2=0;     
    %% compute snell_a2 using total reflection law or snellslaw: 
    if abs(sin(snell_a1/180*pi)*snell_ri1/snell_ri2)>1 
		%%total reflection
		snell_a2=180-snell_a1;
        totalreflection=1;
		%fprintf('\t snellslaw - total reflection, snell_a1: %.2f, snell_a2: %.2f\n\t\t',snell_a1,snell_a2)
    else
		%% regular positive refraction
        snell_a2 = 180/pi*asin(sin(snell_a1/180*pi)*snell_ri1/snell_ri2);
       % snell_a2 = snell_a2+180;
    end %if  
    
    %% determine snell_d2;
    %! does this work for all situation e.g. downward rays, negative angles etc?
    if (snell_a2<-90 || snell_a2>90) && not(totalreflection)
        fprintf('snell_a2 out of bound %.3f\n',snell_a2);
    end
    
    if  surfaceangle<snell_d1
        snell_d2=surfaceangle+270-180-snell_a2;
    else
        snell_d2=surfaceangle+90-180-snell_a2;
    end
    
    if down
        down=1;
        snell_d1=snell_d1+180;
        snell_d2=snell_d2+180;
    end
    
    if snell_d2>360
	   %fprintf('\t snell_d2 > 360 deg, snell_d2: %.2f\n',snell_d2)
	   snell_d2=snell_d2-360;
    end
    if snell_d2<0
	   %fprintf('\t snell_d2 < 0 deg, snell_d2: %.2f\n',snell_d2)
	   snell_d2=snell_d2+360;
    end
   
%fprintf('oldir %.2f newdir %.2f sa %.2f sna1 %.2f sna2 %.2f old_ri %.5f new_ri %.5f  \n',snell_d1, snell_d2, surfaceangle,snell_a1,snell_a2,snell_ri1,snell_ri2)
%fprintf('oldir %.2f newdir %.2f sa %.2f sna1 %.2f sna2 %.2f old_ri %.5f new_ri %.5f ri_base %.5f \n',snell_d1, snell_d2, surfaceangle,snell_a1,snell_a2,snell_ri1,snell_ri2,ri_snellbase)

end %% function snellslaw

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function surfaceangle=isolinefinding(domain,snell_po_base,ss)

%% 	find electron iso surface 
%%
%% input
%%
%% output
%%
%% remarks
%! experiment with different/adaptive circleradius and number of perimeter angles
%! do something about 1) totally flat surfaces 2) strange surfaces (,e.g. interpolation)
%% internal parameters
    circleradius=ss/7;
    numberofcircledivisions=360;

    %% draw circle with radius ss around po_snellbase and determine ri
    perimeter=zeros(numberofcircledivisions,4);
    perimeter(:,1)=1:360/numberofcircledivisions:numberofcircledivisions;    %angles
    
    ri_snellbase=interpolation(domain,snell_po_base,domain.nova-1);

    for a=1:size(perimeter,1)
	   perimeter(a,2)=snell_po_base(1)+cos(perimeter(a,1)/180*pi)*circleradius; %x
	   perimeter(a,3)=snell_po_base(2)+sin(perimeter(a,1)/180*pi)*circleradius; %y
	   perimeter(a,4)=interpolation(domain,perimeter(a,2:3),domain.nova-1); %ri

	   %this is sometimes discontinous, likely poor interpolation on a poor mesh effect
	   %!check it        
    end
    
    %!check uniformity amd continuity for ri(over circle)            
    %for debugging:
	  % plot(perimeter(:,1),perimeter(:,4));
	  % pause
        
    %%find on perimeter the two points that are closest to ri i.e. find
    %%indizes of perimeter, this defines the ri-isosurface
    a1=1;a2=360;    
    diffri1 =abs(perimeter(1,4)-ri_snellbase);
    for a=1:1:179
		if abs(perimeter(a,4)-ri_snellbase)<diffri1
			diffri1=abs(perimeter(a,4)-ri_snellbase);
			a1=a;
		end
    end
 
    diffri2 =abs(perimeter(360,4)-ri_snellbase);
    for a=359:-1:180
		if abs(perimeter(a,4)-ri_snellbase)<diffri2
			diffri2=abs(perimeter(a,4)-ri_snellbase);
			a2=a;
		 end
    end
    surfaceangle=a1;

%{
    %as of 25th sept 2017 all warnings/errors have been verified not to impact, therefore the
    %following check is disabled.
    %% check flatness of surface, warn or terminate depending on level non-flatness
    %% if abs(a2-a1)>356 then the determination of the isosurface fails
    %% likely because ri is constant over this circle 
%    if (abs(a2-a1)<170 || abs(a2-a1)>190) && abs(a2-a1)<356 
        p = polyfit(perimeter(:,1),perimeter(:,4),4)
        %[p,S,mu] = polyfit(x,y,n) also returns mu, which is a two-element vector with centering and scaling values. mu(1) is mean(x), and mu(2) is std(x). Using these values, polyfit centers x at zero and scales it to have unit standard deviation
        fitter=polyval(p,perimeter(:,1));
        wendep1=0.25*p(2)/p(1)+sqrt((0.25*p(2)/p(1))^2-1/6*p(3)/p(1));
        wendep2=0.25*p(2)/p(1)-sqrt((0.25*p(2)/p(1))^2-1/6*p(3)/p(1));
        plot(perimeter(:,1),perimeter(:,4)); %plot ri over angles
        hold on
        t=text(0,ri_snellbase,strcat('a1=',num2str(a1),' a2=',num2str(a2))); %show in figure value of two angles
        plot(perimeter(:,1),fitter);
        warning('snellslaw - ri isosurface not well defined, a1: %.2f, a2: %.2f, a2-a1:%.2f',a1,a2,a2-a1)
        pause
        close;
 %   end %if        
%}
    %{
    if (abs(a2-a1)<140 || abs(a2-a1)>220) && abs(a2-a1)<356 
        %plot(perimeter(:,1),perimeter(:,4)); %plot ri over angles
        %t=text(0,ri_snellbase,strcat('a1=',num2str(a1),' a2=',num2str(a2))); %show in figure value of two angles
		%error('snellslaw - ri isosurface insuffciently defined, a1: %.2f, a2: %.2f, a2-a1:%.2f',a1,a2,a2-a1)
    end %if        
%}
    
end %isoline finding

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function zonenumber=checkifinsidedomain(po,domain)
%%  this function returns the zone and triangle (through the triangle's 3 
%%  vertex coordinates) in which the point is located. Returns zonenumber=0 if the point is outside the domain

%%  input:
%%  double(2)   po          coordinates of point that is checked to be inside domain (which zone and which triangle)
%%  struct      domain      contains all flowfield variables
%%  output:
%%  int         zonenumber  number of zone in which the point po is located
%%  double(2,3) triangle    coordinates of the three vertices of the triangle in which the point po is located
%%  internal parameters
    zonenumber=0;

   %% find closest triangle i.e. its coordinates and indices
   %% convert mesh to vector       
        for z=1:domain.nozones
            if inpolygon(po(1),po(2),  domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  )  
                zonenumber=z;
            end
        end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function value=interpolation(domain,po,index)
%% this function returns an interpolated value

%%  input
%%
%%  output
%%
%%  internal parameters
    value=0;     

    %% find triangle in which the po is inside
    zonenumber=0;
    for z=1:domain.nozones
        %!checkuseofvariables!
        if inpolygon(po(1),po(2),  domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  )  
            zonenumber=z;
        end
    end
    
    if zonenumber
        %% ti = pointLocation(T,QP) returns the enclosing triangles or tetrahedra for each query point in QP. Each row in matrix QP contains the x, y, (and possibly z) coordinates of a query point.
        triangleindex=pointLocation(     domain.( strcat('zone',num2str(zonenumber)) ).DT  ,  po(1),po(2)   );
        pointindices=domain.( strcat('zone',num2str(zonenumber)) ).DT.ConnectivityList( triangleindex,: );

        triangle(:,1)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(1));
        triangle(:,2)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(2));
        triangle(:,3)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(3));

        value= griddata(triangle(1,:),triangle(2,:),triangle(index,:),po(1),po(2),'linear');        
    end
end %% function interpolation



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%%%%source: https://nl.mathworks.com/matlabcentral/fileexchange/38858-contour-plot-for-scattered-data?focused=5249779&tab=function
function [c,h]=tricontour(tri,x,y,z,nv)
%TRICONTOUR Triangular Contour Plot.
% TRICONTOUR(TRI,X,Y,Z,N) draws scalar N contour lines treating the values
% in Z as heights above a plane. TRI,X,Y,and Z define a triangulation where
% the triangles are defined by the M-by-3 face matrix TRI, such as that
% returned by DELAUNAY. Each row of TRI contains indices into the X,Y, and
% Z vertex vectors to define a single triangular face. Contours are
% computed directly from the triangulation rather than interpolating back
% to a cartesian grid using GRIDDATA.
% TRICONTOUR(TRI,X,Y,Z,V) draws length(V) contour lines at the values
% specified in vector V.
% TRICONTOUR(TRI,X,Y,Z,[v v]) draws a single contour line at the level v.
%
% [C,H] = TRICONTOUR(...) returns contour matrix C as described in CONTOURC
% and a vector of handles H to the created patch objects.
% H can be used to set patch properties.
% CLABEL(C) or CLABEL(C,H) labels the contour levels.
%
% Example:
%           x=linspace(-3,3,39);
%           y=linspace(-2.5,2.5,49);
%           [xx,yy]=meshgrid(x,y);
%           zz=peaks(xx,yy);
%           v=-3:1:5; % contour levels
%           subplot(1,2,1)
%           [C,h]=contour(xx,yy,zz,v);   % standard contour for comparison
%           clabel(C)
%           title Contour
% 
%           idx=randperm(numel(zz));     % grab some scattered indices
%           n=idx(1:ceil(numel(zz)/2))'; % one half of them
%           x=xx(n);                     % get scattered data
%           y=yy(n);
%           z=zz(n);
%           tri=delaunay(x,y);           % triangulate scattered data
%           subplot(1,2,2)
%           [C,h]=tricontour(tri,x,y,z,v);
%           clabel(C,h)
%           title TriContour
%
% view(3) displays the contour in 3-D.
%
% See also DELAUNAY, CONTOUR, TRIMESH, TRISURF, TRIPLOT, PATCH.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2006-05-07, 2006-05-16, 2006-07-25

    if nargin<5
        error('Not Enough Input Arguments.')
    end
    x=x(:);	% convert input data into column vectors
    y=y(:);
    z=z(:);
    xlen=length(x);
    if ~isequal(xlen,length(y),length(z))
       error('X, Y, and Z Must Have the Same Number of Elements.')
    end
    if size(tri,2)~=3 || any(tri(:)<0) || any(tri(:)>xlen)
       error('TRI Must Be a Valid Triangulation of the Data in X, Y, Z.')
    end

    zs=z(tri);
    zmax=max(max(zs));              % find max and min in z data that is in tri
    zmin=min(min(zs));

    if length(nv)==1                                 % nv is number of contours
       zlev=linspace(zmax,zmin,nv+2);
    elseif length(nv)==2 && nv(1)==nv(2)              % nv is one contour level
       zlev=nv(1);
    else                                       % nv is vector of contour levels
       zlev=sort(nv,'descend');
    end
    zlev(zlev>=zmax | zlev<=zmin)=[];  % eliminate contours outside data limits
    nlev=length(zlev);

    if nlev==0
       error('No Contours to Plot. Chosen Contours Outside Limits of Data.')
    end

    % precondition the input data
    [zs,zidx]=sort(zs,2);         % sort vertices by z value ascending
    for k=1:size(zs,1)            % shuffle triangles to match
       tri(k,:)=tri(k,zidx(k,:));
    end

    hax=newplot;                  % create new axis if needed
    h=[];                         % patch handle storage
    C=zeros(2,0);                 % Clabel data storage
    cs=[2 1];                     % column swap vector cs(1)=2, cs(2)=1;

    % Main Loop ---------------------------------------------------------------
    for v=1:nlev                  % one contour level at a time
       zc=zlev(v);                % chosen level
       above=zs>=zc;              % true for vertices above given contour
       numabove=sum(above,2);     % number of triangle vertices above contour
       tri1=tri(numabove==1,:);   % triangles with one vertex above contour
       tri2=tri(numabove==2,:);   % triangles with two vertices above contour
       n1=size(tri1,1);           % number with one vertex above
       n2=size(tri2,1);           % number with two vertices above

       edge=[tri1(:,[1 3])        % first column is indices below contour level
             tri1(:,[2 3])        % second column is indices above contour level
             tri2(:,[1 2])
             tri2(:,[1 3])];
       if n1==0                   % assign edges to triangle number
          n=[1:n2 1:n2]';
       elseif n2==0
          n=[1:n1 1:n1]';
       else
          n=[1:n1 1:n1 n1+(1:n2) n1+(1:n2)]';
       end

       [edge,idx]=sortrows(edge);    % put shared edges next to each other
       n=n(idx);                     % shuffle triangle numbers to match

       idx=all(diff(edge)==0,2);     % find shared edges
       idx=[idx;false]|[false;idx];  % True for all shared edges

       % eliminate redundant edges, two triangles per interior edge
       edgeh=edge(~idx,:);           % hull edges
       nh=n(~idx);                   % hull triangle numbers
       if ~isempty(nh)
          nh(end,2)=0;               % zero second column for hull edges
       end
       edges=edge(idx,:);            % shared edges
       edges=edges(1:2:end-1,:);     % take only unique edges
       ns=n(idx);                    % interior triangle numbers
       ns=[ns(1:2:end) ns(2:2:end)]; % second column is second triangle
       edge=[edgeh;edges];           % unique edges
       nn=[nh;ns];                   % two columns of triangle numbers
       ne=size(edge,1);              % number of edges

       flag=true(ne,2);              % true for each unused edge per triangle
       tmp=zeros(ne+1,1);            % contour data temporary storage

       xe=x(edge);                   % x values at vertices of edges
       ye=y(edge);                   % y values at  vertices of edges
       ze=z(edge);                   % z data at  vertices of edges

       alpha=(zc-ze(:,1))./(ze(:,2)-ze(:,1)); % interpolate all edges
       xc=alpha.*(xe(:,2)-xe(:,1)) + xe(:,1); % x values on this contour
       yc=alpha.*(ye(:,2)-ye(:,1)) + ye(:,1); % y values on this contour

       while any(flag)	% while there are still unused edges -----------------

          xtmp=tmp;
          ytmp=tmp;
          [ir,ic]=find(flag,1);            % find next unused edge
          flag(ir,ic)=false;               % mark this edge used

          k=1;                             % first data point in subcontour
          xtmp(k)=xc(ir);                  % store data from this edge
          ytmp(k)=yc(ir);

          while true     % complete this subcontour ---------------------------

             [ir,ic]=find(flag&nn(ir,ic)==nn,1);% find other edge of triangle
             flag(ir,ic)=false;            % mark this edge used
             k=k+1;
             xtmp(k)=xc(ir);               % store data from this edge
             ytmp(k)=yc(ir);

             ic=cs(ic);                    % other triangle that shares edge

             if nn(ir,ic)==0               % reached hull, subcontour complete
                k=k+1;
                xtmp(k)=nan;               % don't let subcontour close
                ytmp(k)=nan;
                break
             elseif ~flag(ir,ic)           % complete closed subcontour
                break
             else                          % more points remain on subcontour
                flag(ir,ic)=false;         % mark this edge used
             end
          end % while true ----------------------------------------------------
          xtmp(k+1:end)=[];                % throw away unused storage
          ytmp(k+1:end)=[];                % xtmp,ytmp contain subcontour

          if nargout<2                     % plot the subcontour
             patch('XData',xtmp,'YData',ytmp,'CData',repmat(zc,k,1),...
                   'Parent',hax,'FaceColor','none','EdgeColor','flat',...
                   'UserData',zc)
             C=horzcat(C,[zc xtmp';k ytmp']); % contour label data
          else                             % plot subcontour and create output
             h=[h;patch('XData',xtmp,'YData',ytmp,'CData',repmat(zc,k,1),...
             'Parent',hax,'FaceColor','none','EdgeColor','flat',...
             'UserData',zc)]; %#ok
             C=horzcat(C,[zc xtmp';k ytmp']); % contour label data
          end
       end % while any(flag) --------------------------------------------------
    end % for v=1:nlev
    if nargout
       c=C;
    end
end


%{
-----------------------------------------------------------------------
Number of species: 
11
-----------------------------------------------------------------------
Species symbols: 
em
N
O
N2
NO 
O2
N2p
NOp
Np
O2p
Op
-----------------------------------------------------------------------
STOP

%}




%pltformat 1
%{ 
TITLE     = "| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |"
VARIABLES = "x0"
"x1"
"numb_e"
"numb_CO2"
"numb_N2"
"numb_C"
"numn_N"
"numb_O"
"numb_O2"
"numb_CO"
"numb_NO"
"numb_CO+"
"numb_NO+"
"numb_C+"
"numb_O+"
"numb_O2+"
"u"
"v"
"T_LARS"
"rho"
"H"
"M"
"p"
ZONE T="ExoMars_t_73s"
 STRANDID=1, SOLUTIONTIME=0
 Nodes=77900, Elements=77244, ZONETYPE=FEQuadrilateral
 DATAPACKING=POINT
 DT=(SINGLE SINGLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE SINGLE SINGLE DOUBLE SINGLE SINGLE SINGLE SINGLE )
 4.402699172E-01 1.199999809E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 1.600104151E-03 6.285274500E+06 0.000000000E+00 4.299663696E+02
 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 4.123243690E-02 6.549364000E+06 0.000000000E+00 1.130732910E+04
 1.578437388E-01 1.852897167E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 4.516700195E+03 2.517777653E-14 0.000000000E+00 5.820000079E-04 1.031773200E+07 2.130347061E+01 1.924182701E+01
 -1.694117188E-01 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 4.516701172E+03 -6.373672513E-04 0.000000000E+00 5.819887738E-04 1.031773700E+07 2.130354309E+01 1.924132156E+01
 
%}

%pltformat 2
%{

TITLE     = "| | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |"
VARIABLES = "x0" "x1" "numb_e" "numb_CO2" "numb_N2" "numb_C" "numn_N" "numb_O" "numb_O2" "numb_CO" "numb_NO" "numb_CO+" "numb_NO+" "numb_C+" "numb_O+" "numb_O2+" "u" "v" "T_LARS" "rho" "H" "M" "p"
ZONE T="ExoMars_t_73s"  N=77900, E=77244, F=FEPOINT, ET=QUADRILATERAL, AUXDATA CPU="0", AUXDATA TRS="InnerCells", AUXDATA Filename="Star.plt", AUXDATA ElementType="Quad", AUXDATA Iter="1", AUXDATA PhysTime="0"
 4.402699172E-01 1.199999809E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 1.600104151E-03 6.285274500E+06 0.000000000E+00 4.299663696E+02
 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 4.123243690E-02 6.549364000E+06 0.000000000E+00 1.130732910E+04
 1.578437388E-01 1.852897167E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 4.516700195E+03 2.517777653E-14 0.000000000E+00 5.820000079E-04 1.031773200E+07 2.130347061E+01 1.924182701E+01
 -1.694117188E-01 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 4.516701172E+03 -6.373672513E-04 0.000000000E+00 5.819887738E-04 1.031773700E+07 2.130354309E+01 1.924132156E+01
 1.321143270E+00 9.481133893E-03 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 0.000000000E+00 9.269155562E-04 6.680586000E+06 0.000000000E+00 2.553797913E+02

%}



