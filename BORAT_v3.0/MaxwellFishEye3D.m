clc; clear all; close all;

runFolder = 'Solutions';


%domain size
x=linspace(-3,3,100);
y=linspace(-3,3,100);
z=linspace(-3,3,100);

%radius of the lens
R=2;

[X,Y,Z]=meshgrid(x,y,z);

n=zeros(length(x),length(y),length(z));
nx=zeros(length(x),length(y),length(z));
ny=zeros(length(x),length(y),length(z));
nz=zeros(length(x),length(y),length(z));

for k=1:length(y)
    for j=1:length(x)
        for l=1:length(z)
            Rk=sqrt((X(k,j,l).^2+Y(k,j,l).^2+Z(k,j,l)^2));
            
            if Rk>R
                n(k,j,l)=1;
                nx(k,j,l)=0;
                ny(k,j,l)=0;
                nz(k,j,l)=0;
            else
                n(k,j,l)=2/(1+(Rk/R)^2);
                nx(k,j,l)=(-4*(X(k,j,l)/R^2)/((1+(Rk/R)^2)^2));
                ny(k,j,l)=(-4*(Y(k,j,l)/R^2)/((1+(Rk/R)^2)^2));
                nz(k,j,l)=(-4*(Z(k,j,l)/R^2)/((1+(Rk/R)^2)^2));
                
            end
        end
    end
end

grad_n=sqrt(nx.^2+ny.^2+nz.^2);

figure;
isosurface(X,Y,Z,n,1);


%% Eikonal integration

figure;
xlim([-3 3])
ylim([-3 3])
hold on;


alpha1=0;
alpha2=360;

n_alpha=15; %n_point per ring
n_theta=10; %n_ring

Aperture=179;

body_normalZ=90;
body_normalXY=90;

d_alpha=(alpha2-alpha1)/(n_alpha);
d_theta=(Aperture/2)/(n_theta);

for i=1:n_alpha
    
    alpha(i)=  alpha1+d_alpha*(i-1);
    
end

for i=1:n_theta
    
    theta(i)=  90-d_theta*(i);
    
end

itpo=zeros(4,100,n_theta*n_alpha);
nRay=0;
for k=1:length(theta)
    for i=1:length(alpha)
        
        nRay=nRay+1;
        
        x_0=-2;
        y_0=0;
        z_0=0;
        
        alpha_0=alpha(i);
        theta_0=theta(k);
        
        [n_0,nx_0,ny_0,nz_0] = refractive(x_0,y_0,z_0);
        
        
        chi_x_0=n_0*(cos(theta_0*pi/180))*(cos(alpha_0*pi/180));
        chi_y_0=n_0*(cos(theta_0*pi/180))*(sin(alpha_0*pi/180));
        chi_z_0=n_0*sin(theta_0*pi/180);
        
        chi2=roty(body_normalZ)*[chi_x_0 chi_y_0 chi_z_0]';
        
        y0=[x_0,y_0,z_0,chi2(1),chi2(2),chi2(3)];
        
        tspan = [0 100];
        
        % [t,y_solution] = ode45(@eikonal,tspan,y0);
        tstart = 0;
        tfinal = 10;
        
        refine = 4;
        
        options = odeset('Events',@cutoff,...
            'OutputSel',1,'Refine',refine);
        
        options = odeset(options,'MaxStep',0.5);
        % options =odeset(options, 'OutputFcn', @odeplot);
        
        tout = tstart;
        yout = y0;
        teout = [];
        yeout = [];
        ieout = [];
        for i = 1:1
            
            
            [t,y,te,ye,ie] = ode45(@eikonal,[tstart tfinal],y0,options);
            
            
            if ~ishold
                hold on
            end
            % Accumulate output.
            nt = length(t);
            tout = [tout; t(2:nt)];
            yout = [yout; y(2:nt,:)];
            teout = [teout; te];    % Events at tstart are never reported.
            yeout = [yeout; ye];
            ieout = [ieout; ie];
            
            % A good guess of a valid first time step is the length of
            % the last valid time step, so use it for faster computation.
            options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
                'MaxStep',t(nt)-t(1));
            
            
            tstart = t(nt);
        end
        
        
        
        plot3(yout(:,1),yout(:,2), yout(:,3),'-b','Linewidth',2);
        
        itpo(1,1:length(yout(:,1)),nRay)= yout(:,1); % x
        itpo(2,1:length(yout(:,1)),nRay)= yout(:,2); % y
        itpo(3,1:length(yout(:,1)),nRay)= yout(:,3); % y
        
        for  j=1:length(yout(:,1))
            itpo(4,j,nRay)= sqrt(yout(j,4)^2+yout(j,5)^2+yout(j,6)^2); % refractive along path (not integrated)
        end
        
    end
end
hold off

% %% Export solutions tecplot
% 
% %path for tecplot export files
% 
% tecplot_home = '/Applications/Tecplot 360 EX 2020 R2/Tecplot 360 EX 2020 R2.app/Contents';
% tecio_path = strcat(tecplot_home, '/Frameworks/libtecio.dylib');
% 
% tecio_header_path = '/TecIO/TECIO.h'; %Uses simplified TECIO.h included in package.
% 
% output_fname = 'Maxwell3D_domain.szplt';
% output_fname_rays = 'Maxwell_Rays_3D.szplt';
% 
% I = size(n,1); J = size(n,2); K = size(n,3);
% total_points = I*J*K;
% 
% Xe = reshape(X,total_points,1);
% Ye = reshape(Y,total_points,1);
% Ze = reshape(Z,total_points,1);
% Ne = reshape(n,total_points,1);
% GradNe = reshape(grad_n,total_points,1);
% vars = [Xe,Ye,Ze,Ne,GradNe];
% 
% if ~libisloaded('tecio')
%     [notfound, warnings] = loadlibrary(tecio_path, tecio_header_path,...
%         'alias', 'tecio');
% end
% 
% 
% dataset_title = 'Maxwell Fisheye';
% var_names = 'X, Y, Z, n, gradN';
% file_format = 1;  % 0 - .plt;  1 - .szplt
% file_type   = 0;  % 0 - grid & solution; 1 - grid only; 2 - solution only
% data_type   = 2;  % 1 - single; 2 - double; ...
% 
% [isok,~,~,~,~,filehandle] = calllib('tecio', 'tecFileWriterOpen', ...
%     output_fname, dataset_title, var_names, ...
%     file_format, file_type, data_type, [],[] );
% 
% zname = dataset_title;
% z_idx = 1;
% [isok,~,~,~,~,~,~,z_idx] = calllib('tecio', 'tecZoneCreateIJK', ...
%     filehandle, zname, I, J, K, [], [], [],[],0,0,0,z_idx);
% 
% 
% for v = 1:size(vars,2)
%     isok = calllib('tecio', 'tecZoneVarWriteDoubleValues', ...
%         filehandle, z_idx, v, 1, total_points, vars(:,v));
% end
% 
% 
% calllib('tecio', 'tecFileWriterClose', filehandle);
% 
% if libisloaded('tecio')
%     unloadlibrary('tecio')
% end
% 
% % Export rays Eikonal
% 
% 
% dataset_title = 'Ray solution ';
% var_names = 'X, Y, Z, n';
% file_format = 1;  % 0 - .plt;  1 - .szplt
% file_type   = 0;  % 0 - grid & solution; 1 - grid only; 2 - solution only
% data_type   = 2;  % 1 - single; 2 - double; ...
% 
% total_zones=size(itpo(1,1,:),3);
% 
% if ~libisloaded('tecio')
%     [notfound, warnings] = loadlibrary(tecio_path, tecio_header_path,...
%         'alias', 'tecio');
% end
% 
% zone_title = 'Ray n ';
% 
% [isok,~,~,~,~,filehandle] = calllib('tecio', 'tecFileWriterOpen', ...
%     output_fname_rays, dataset_title, var_names, ...
%     file_format, file_type, data_type, [],[] );
% 
% wzone=0;
% for nzone=1:1:(total_zones)
%     wzone=wzone+1;
%     total_points=nnz(itpo(1,:,nzone));
%     I=total_points;
%     J=1;
%     K=1;
%     zname = [ zone_title num2str(nzone,'%02d')];
%     
%     [isok,~,~,~,~,~,~,n_zone] = calllib('tecio', 'tecZoneCreateIJK', ...
%         filehandle, zname, I, J, K, [], [], [],[],0,0,0,wzone);
%     
%     vars= [itpo(1,1:total_points,nzone)' itpo(2,1:total_points,nzone)' itpo(3,1:total_points,nzone)' itpo(4,1:total_points,nzone)' ];
%     
%     for v = 1:size(vars,2)
%         isok = calllib('tecio', 'tecZoneVarWriteDoubleValues', ...
%             filehandle, wzone, v, 1, total_points, vars(:,v));
%     end
%     
% end %nzone
% 
% calllib('tecio', 'tecFileWriterClose', filehandle);
% 
% if libisloaded('tecio')
%     unloadlibrary('tecio')
% end
% 
% command=strcat('mv *.szplt',32, runFolder); %32 is the space code
% unix(command);

%% Functions for eikonal solver

function [value,isterminal,direction] = cutoff(t,y)
    
    
    [n, nx, ny, nz] = refractive(y(1), y(2), y(3));
    flag=1;
    if n<0.1
        flag=0;
    end
    
    if abs(y(1))>2
        flag=0;
    end
    
    if abs(y(2))>2
        flag=0;
    end
    
    if abs(y(3))>2
        flag=0;
    end
    
    value = flag;
    isterminal = 1;            % Stop at local minimum
    direction = -1;            % [local minimum, local maximum]
end   % End nested function events



function [n,nx,ny,nz] = refractive(x,y,z)
    
    
    R=2;
    Rk=sqrt(x^2+y^2+z^2);
    
    if Rk>R
        n=1;
        ny=0;
        nx=0;
        nz=0;
    else
        n=2/(1+(Rk/R)^2);
        nx=(-4*x/R^2)/((1+(Rk/R)^2)^2);
        ny=(-4*y/R^2)/((1+(Rk/R)^2)^2);
        nz=(-4*z/R^2)/((1+(Rk/R)^2)^2);
        
    end
end



function solution = eikonal(t,y)
    
    
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
