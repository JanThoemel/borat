clc; clear all; close all;

n_alpha=15; % N longitudinal rays

flowfieldfilename='3Dmeshes/Paraboloid3D.dat';

wallfilename='3Dmeshes/ParaboloidWall3D.dat';

%% Read domain file

z=0;
domain.nova=3; %just 3D points
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




figure;
patch('Vertices',domain.( strcat('zone',num2str(z)) ).BoundaryPoints,...
    'Faces',domain.( strcat('zone',num2str(z)) ).BoundaryFaces,...
    'FaceColor','w','FaceAlpha',1,'EdgeColor','b');
axis equal
view(51,24)

%% Read wall file

z=0;
domain.nova=3; %just 3D points
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
    
    domain.( strcat('zone',num2str(z)) ).Wall.variables=zeros(  domain.nova  ,  domain.( strcat('zone',num2str(z)) ).Wall.N);
    
    
    
    for i=1:domain.( strcat('zone',num2str(z)) ).Wall.N
        line=fgetl(id);
        domain.( strcat('zone',num2str(z)) ).Wall.variables(1:domain.nova,i)=sscanf(line,'  %g');
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






figure;
patch('Vertices',domain.( strcat('zone',num2str(z)) ).Wall.variables(1:3,:)',...
    'Faces',domain.( strcat('zone',num2str(z)) ).Wall.NodeMap,...
    'FaceColor','w','FaceAlpha',1,'EdgeColor','b');
axis equal
view(-51,24)



%% Eikonal integration

figure;hold on;
patch('Vertices',domain.( strcat('zone',num2str(z)) ).Wall.variables(1:3,:)',...
    'Faces',domain.( strcat('zone',num2str(z)) ).Wall.NodeMap,...
    'FaceColor','b','FaceAlpha',.7,'EdgeColor','k');
view(0,0)

alpha1=0;
alpha2=360;

n_theta=1; %

d_alpha=(alpha2-alpha1)/(n_alpha);


for i=1:n_alpha
    
    alpha(i)=  alpha1+d_alpha*(i-1);
    
end


theta=0;

nResemple=50;

itpo=zeros(4,100,n_theta*n_alpha);
nRay=0;

for k=1:length(theta)
    for i=1:length(alpha)
        
        nRay=nRay+1;
        
        x_0=0;
        y_0=0;
        z_0=0.025;
        
        alpha_0=alpha(i);
        theta_0=theta(k);
        
        [n_0,nx_0,ny_0,nz_0] = refractive(x_0,y_0,z_0);
        
        
        chi_x_0=n_0*(cos(theta_0*pi/180))*(cos(alpha_0*pi/180));
        chi_y_0=n_0*(cos(theta_0*pi/180))*(sin(alpha_0*pi/180));
        chi_z_0=n_0*sin(theta_0*pi/180);
        
        y0=[x_0,y_0,z_0,chi_x_0,chi_y_0,chi_z_0];
        tspan = [0 100];
        
        tstart = 0;
        tfinal = 10;
        
        refine = 4;
        
        options = odeset('Events',@(t,y) cutoff(t,y,domain));
        
        
        options = odeset(options,'Refine',refine); %this doesn'' make it work
        
        options = odeset(options,'MaxStep',0.001);
        
        tout = tstart;
        yout = y0;
        teout = [];
        yeout = [];
        ieout = [];
        
        y_final=[y0(1) y0(2) y0(3)];
        
        flag_out=checkifinsidedomain3D(y_final, domain);
        
        while flag_out
            
            solution=ode45(@(t,y) eikonal(t,y,domain),[tstart tfinal],y0,options);
            
            
            t_sol=linspace(tstart,solution.x(end),nResemple);
            
            ray_solution=deval(solution,t_sol);
            
            % Accumulate output.
            nt = length(solution.y);
            tout = [tout; solution.x(1,2:nt)'];
            yout = [yout; solution.y(:,2:nt)'];
            teout = [teout; solution.xe];    % Events at tstart are never reported.
            yeout = [yeout; solution.ye'];
            ieout = [ieout; solution.ie];
            
            
            flag_boundary= checkBoundaryReflection3D([yout(end,1) yout(end,2) yout(end,3)],domain);
            
            if flag_boundary==1
                
                P1=[yout(end-1,1) yout(end-1,2) yout(end-1,3)];
                P2=[yout(end,1) yout(end,2) yout(end,3)]; %this is the new starting point
                u=P2-P1;
                
                
                IDboundary=nearestNeighbor(domain.(strcat('zone',num2str(z))).Wall.BoundaryDT,P2);
                N2=vertexNormal(domain.(strcat('zone',num2str(z))).Wall.BoundaryDT,IDboundary);
                
                IDpointCenter=dsearchn(domain.( strcat('zone',num2str(z)) ).Wall.CenterPoint, P2);
                N=faceNormal(domain.(strcat('zone',num2str(z))).Wall.BoundaryDT,IDpointCenter);
                
                if N(1)~=N2(1) && N(2)~=N2(2) && N(3)~=N2(3)
                    N;
                    N2;
                end
                
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
                
                
                [n_0,nx_0,ny_0,nz_0] = refractive(yout(end,1),yout(end,2), yout(end,3));
                
                y0(4)=n_0*cos(theta0*pi/180)*cos(alpha0*pi/180);
                y0(5)=n_0*cos(theta0*pi/180)*sin(alpha0*pi/180);
                y0(6)=n_0*sin(theta0*pi/180);
                
                y0(1) = yout(end-2,1);
                y0(2) = yout(end-1,2);
                y0(3) = yout(end-1,3);
                
                options = odeset(options,'InitialStep',tout(nt)-tout(nt-refine));
                
                tstart = tout(nt);
                
                plot3(yout(:,1),yout(:,2), yout(:,3),'-r','Linewidth',1);
                hold on;
                
            else
                
                flag_out=0;
            end
            
            options = odeset(options,'InitialStep',tout(nt)-tout(nt-refine),...
                'MaxStep',tout(nt)-tout(1));
            
            tstart = tout(nt);
            plot3(yout(:,1),yout(:,2), yout(:,3),'-r','Linewidth',2);
            hold on;
            
        end
        
        itpo(1,1:length(yout(:,1)),nRay)= yout(:,1); % x
        itpo(2,1:length(yout(:,1)),nRay)= yout(:,2); % y
        itpo(3,1:length(yout(:,1)),nRay)= yout(:,3); % y
        
        for  j=1:length(yout(:,1))
            itpo(4,j,nRay)= sqrt(yout(j,4)^2+yout(j,5)^2+yout(j,6)^2); % refractive along path (not integrated)
        end
        
    end
end
axis equal
hold off



%% Functions for eikonal solver

function [value,isterminal,direction] = cutoff(t,y,domain)
    
    flag_out=checkifinsidedomain3D([y(1) y(2) y(3)],domain);
    
    if flag_out~=1
        
        value = flag_out;
        isterminal = 1;            % Stop at local minimum
        direction = 0;
        
    else
        
        [n, nx, ny, nz] = refractive(y(1), y(2), y(3));
        flag=1;
        if n<0.1
            flag=0;
        end
        
        if abs(y(1))>abs(2)
            flag=0;
        end
        
        if abs(y(2))>(2)
            flag=0;
        end
        
        if abs(y(3))>(1)
            flag=0;
        end
        
        value = flag;
        isterminal = 1;            % Stop at local minimum
        direction = -1;            % [local minimum, local maximum]
    end
    
end   % End nested function events



function [n,nx,ny,nz] = refractive(x,y,z)
    
    n=1;
    ny=0;
    nx=0;
    nz=0;
    
end



function solution = eikonal(t,y,domain)
    
    
    xi=y(1);
    yi=y(2);
    zi=y(3);
    ui=y(4);
    vi=y(5);
    wi=y(6);
    
    
    [n, nx, ny, nz] = refractive(xi, yi, zi);
    
    
    yp = zeros(6,1);
    
    yp(1) = ui / (n^2) ;
    yp(2) = vi / (n^2) ;
    yp(3) = wi / (n^2) ;
    
    yp(4) = nx/n ;
    yp(5) = ny/n ;
    yp(6) = nz/n ;
    
    solution = [yp(1); yp(2); yp(3); yp(4); yp(5); yp(6)];
    
end