clc; clear all; close all;

runFolder = 'Solutions';


%domain size
x=linspace(-3,3,100);
y=linspace(-3,3,100);

%radius of the lens
R=2;


[X,Y]=meshgrid(x,y);

n=zeros(length(y),length(x));
nx=zeros(length(y),length(x));
ny=zeros(length(y),length(x));




for k=1:length(y)
    for j=1:length(x)
        Rk=sqrt((X(k,j).^2+Y(k,j).^2));
        
        if Rk>R
            n(k,j)=1;
            ny(k,j)=0;
            nx(k,j)=0;
        else
            n(k,j)=2/(1+(Rk/R)^2);
            nx(k,j)=(-4*X(k,j)/R^2)/((1+(Rk/R)^2)^2);
            ny(k,j)=(-4*Y(k,j)/R^2)/((1+(Rk/R)^2)^2);
        end
    end
end

grad_n=sqrt(nx.^2+ny.^2);

figure;
[C,h]=contourf(X,Y,n,10);
colormap jet
colorbar
caxis([1 1.5])
title('refractive index')


%%%%%%%%%%%%%%%


%% Eikonal integration

figure;
%surf(X,Y,n);
[C,h]=contourf(X,Y,n,100);
set(h,'Linecolor','none');
xlim([-3 3])
ylim([-3 3])
hold on;
axis equal

n_ray=30;

alpha_0_vec=linspace(0,360,n_ray);

itpo=zeros(3,100,n_ray);

for k=1:n_ray
    
    x_0=-2;
    y_0=0;
    
    alpha_0=alpha_0_vec(k);
    
    [n_0,nx_0,ny_0] = refractive(x_0,y_0);
    
    chi_x_0=n_0*cos(alpha_0*pi/180);
    chi_y_0=n_0*sin(alpha_0*pi/180);
    
    
    y0=[x_0,y_0,chi_x_0,chi_y_0];
    tspan = [0 10];
    
    
    tstart = 0;
    tfinal = 10;
    
    y0 = [x_0,y_0,chi_x_0,chi_y_0];
    
    refine = 4;
    
    options = odeset('Events',@cutoff,...
        'OutputSel',1,'Refine',refine);
    
    options = odeset(options,'MaxStep',0.5);
    
    
    tout = tstart;
    yout = y0;
    teout = [];
    yeout = [];
    ieout = [];
    for i = 1:1
        % Solve until the first terminal event.
        
        
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
        
        % Set the new initial conditions to the step back mirrored.
        y0(1) = y(nt,1);
        y0(2) = y(nt,2);
        y0(3) = y(nt,3);
        y0(4) = -y(nt,4);
        
        % A good guess of a valid first time step is the length of
        % the last valid time step, so use it for faster computation.
        options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
            'MaxStep',t(nt)-t(1));
        
        
        tstart = t(nt);
    end
    
    
    
    plot(yout(:,1),yout(:,2),'k','Linewidth',2);
    
    itpo(1,1:length(yout(:,1)),k)= yout(:,1); % x
    itpo(2,1:length(yout(:,1)),k)= yout(:,2); % y
    
    for  i=1:length(yout(:,1))
        itpo(3,i,k)= sqrt(yout(i,3)^2+yout(i,4)^2); % refractive along path (not integrated)
    end
    
end

hold off
axis equal
% plot(yeout(1),yeout(2),'o','Linewidth',2);

% %% Export solutions tecplot
% 
% %path for tecplot export files
% tecplot_home = '/Applications/Tecplot 360 EX 2020 R2/Tecplot 360 EX 2020 R2.app/Contents';
% tecio_path = strcat(tecplot_home, '/Frameworks/libtecio.dylib');
% tecio_header_path = '/TecIO/TECIO.h'; %Uses simplified TECIO.h included in package.
% 
% output_fname = 'MaxwellFishEye2D_domain.szplt';
% output_fname_rays = 'MaxwellFishEye2D_Rays_2D.szplt';
% 
% 
% I = size(n,1); J = size(n,2); K = 1;
% total_points = I*J*K;
% 
% Xe = reshape(X,total_points,1);
% Ye = reshape(Y,total_points,1);
% Ne = reshape(n,total_points,1);
% GradNe = reshape(grad_n,total_points,1);
% vars = [Xe,Ye,Ne,GradNe];
% 
% if ~libisloaded('tecio')
%     [notfound, warnings] = loadlibrary(tecio_path, tecio_header_path,...
%         'alias', 'tecio');
% end
% 
% 
% dataset_title = 'Maxwell Fisheye';
% var_names = 'X, Y, n, GradNe';
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
% %export ray solution
% 
% dataset_title = 'Ray solution ';
% var_names = 'X, Y, n';
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
%     vars= [itpo(1,1:total_points,nzone)' itpo(2,1:total_points,nzone)' itpo(3,1:total_points,nzone)' ];
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
    
    [n, nx, ny] = refractive(y(1), y(2));
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
    
    
    value = flag;
    isterminal = 1;            % Stop at local minimum
    direction = -1;            % [local minimum, local maximum]
end   % End nested function events


function [n,nx,ny] = refractive(x,y)
    
    
    R=2;
    Rk=sqrt(x^2+y^2);
    
    if Rk>R
        n=1;
        ny=0;
        nx=0;
    else
        n=2/(1+(Rk/R)^2);
        nx=(-4*x/R^2)/((1+(Rk/R)^2)^2);
        ny=(-4*y/R^2)/((1+(Rk/R)^2)^2);
    end
end



function solution = eikonal(t,y)
    
    
    xi=y(1);
    yi=y(2);
    zi=y(3);
    wi=y(4);
    
    
    [n, nx, ny] = refractive(xi, yi);
    
    
    yp = zeros(4,1);
    
    yp(1) = zi / (n) ;
    yp(2) = wi / (n) ;
    yp(3) = nx ;
    yp(4) = ny ;
    
    solution = [yp(1); yp(2); yp(3); yp(4)];
    
end
