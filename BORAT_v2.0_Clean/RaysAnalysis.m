clc; clear all; close all;

set(groot,'DefaultFigureColormap',gray)


load('Solutions/ExoMars/ExoMars_38s.mat');

runFolder='ExoMars_38s_raytracing';
solutionFile='ExoMars_38s_raytracing_solution';

mkdir (runFolder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raytracing parameters


dir1=43+80;					%%ExoMars full case TBC
dir2=43-80;					%%ExoMars full case TBC

pooo1=[1.252+0.001 1.066+0.001];       %% ExoMars 
pooo2=[1.252+0.001 1.066+0.001];       %% ExoMars


maxangles=20;

domain.limitx1=0;
domain.limitx2=1.872;

domain.limity1=-1.7;
domain.limity2=1.7;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Raytracing Solver

%Eikonal
[itdir, itpo,~]=eikonal2D(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, absorptionlimits(2)-10,symmetryline);

draw_raytracing(domain,plotscdir,scdir,scdir_range,collisionmodel,absorptionlimits,pooo1,pooo2,itpo,maxangles,itdir,symmetrylineencounter,writefig,runFolder)


%Snell's law
% [itdir_snell, itpo_snell,symmetrylineencounter]=raytracing(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, absorptionlimits(2)-10,symmetryline,ss);
%draw_raytracing(domain,plotscdir,scdir,scdir_range,collisionmodel,absorptionlimits,pooo1,pooo2,itpo_snell,maxangles,itdir_snell,symmetrylineencounter,writefig,runFolder)


cd (runFolder); 
save(solutionFile);
cd ..;


