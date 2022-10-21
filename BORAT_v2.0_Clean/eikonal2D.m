function [itdir, itpo,symmetrylineencounter]=eikonal2D(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, absorptioncutoff,symmetryline,ss,f)


 %% this function traces rays with eikonal equation
    
    msg=' ';
    symmetrylineencounter=zeros(maxangles,1);
    
    cut_off_limit=0;
    
    %% compute 2D gradient
    
    for z=1:domain.nozones
        
        [domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,:),domain.( strcat('zone',num2str(z)) ).variables(domain.nova+2,:)] = trigradient(domain.( strcat('zone',num2str(z)) ).variables(1,:)',domain.( strcat('zone',num2str(z)) ).variables(2,:)',domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)',domain.( strcat('zone',num2str(z)) ).delaunay);
        
    end
    
    
    %% trace from around wall normal
    if maxangles==1
        anglestep=0;
    else
        anglestep=(dir1-dir2)/(maxangles-1);
        x_step=(pooo1(1)-pooo2(1))/(maxangles-1);
        y_step=(pooo1(2)-pooo2(2))/(maxangles-1);
    end
    
    
    
    %% iteration points: 1: x, 2: y, 3: pathlength, 4: ri, 5: ec
    % here we add the integration of refractive index, that is the
    % optical path lenght
    itpo=zeros(6,maxsteps,maxangles);
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
        %% starting point optical path length: this need to be checked
        itpo(6,1,:)=interpolation(domain,itpo(:,1,1),domain.nova-1);
        %itpo(6,1,:)=0; % this might be an option, but not sure
        
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
            %% starting point optical path length: this need to be checked
            itpo(6,1,:)=interpolation(domain,itpo(:,1,1),domain.nova-1);
            %itpo(6,1,:)=0; % this might be an option, but not sure
            
        end %% for
    end %%if
    
    fprintf('\n  starting point refractive index for all rays %.6f',  itpo(4,1,1));
    fprintf('\n  starting point absorption coefficient for all rays %.6f ',itpo(5,1,1))
    
    %mesh1 and mesh2 contains the point in the wall surface for each
    %zone
    
    
    mesh1=[domain.zone1.variables(1,domain.zone1.bound(:))',domain.zone1.variables(2,domain.zone1.bound(:))'];
    
    if (domain.nozones==2)
        mesh2=[domain.zone2.variables(1,domain.zone2.bound(:))',domain.zone2.variables(2,domain.zone2.bound(:))'];
    end
    
    
    for a=1:maxangles
        
        %% integration initial conditions
        
        x_0=itpo(1,1,a);
        y_0=itpo(2,1,a);
        alpha_0=itdir(1,a);
        
        [n_0,~,~, absorption_0] = refractive_interpolation(x_0,y_0,domain);
        
        chi_x_0=n_0*cos(alpha_0*pi/180);
        chi_y_0=n_0*sin(alpha_0*pi/180);
        
        y0 =[x_0,y_0,chi_x_0,chi_y_0, n_0, absorption_0,0]; %7 terms
        
        %% Ode settings
        
        tstart = 0;
        tfinal =2000;
        
        refine = 4;
        
        options = odeset('Events',@(t,y) cutoff(t,y,domain));
        
        options = odeset(options,'Refine',refine); %this sometimes is problematic
        
        options = odeset(options,'InitialStep',1e-2);
        
        options = odeset(options,'MaxStep',1);
        
        
        tout = tstart;
        yout = y0;
        teout = [];
        yeout = [];
        ieout = [];
        
        
        y_final=[y0(1) y0(2)];
        
        flag_out=checkifinsidedomain([y_final(1) y_final(2)], domain);
        
        msg=sprintf('\n \t a=%d/%d inidir=%.0f/(%.0f - %0.f)',a,maxangles,itdir(1,a),itdir(1,1),itdir(1,maxangles));
        fprintf(msg)
        
        while flag_out
            
            % the tolerance of different solvers gives different results
            % these should be checked better
            
            %[t,y,te,ye,ie] = ode45(@(t,y) eikonal_equation(t,y,domain),[tstart tfinal],y0,options);
            [t,y,te,ye,ie] = ode15s(@(t,y) eikonal_equation(t,y,domain),[tstart tfinal],y0,options);
            %[t,y,te,ye,ie] = ode113(@(t,y) eikonal_equation(t,y,domain),[tstart tfinal],y0,options);
            
            % Accumulate output.
            nt = length(t);
            tout = [tout; t(2:nt)];
            yout = [yout; y(2:nt,:)];
            teout = [teout; te];    % Events at tstart are never reported.
            yeout = [yeout; ye];
            ieout = [ieout; ie];
            
            
            flag_boundary = checkifboundary([yout(end,1) yout(end,2)],domain);
            
            %cutoff check
            flag_cutoff=0;
            [n_end, nx_end, ny_end] = refractive_interpolation(y(end,1), y(end,2), domain);
            
            
            if flag_boundary==1 && n_end>cut_off_limit
                % Set the new initial conditions to the step back mirrored.
                
                fprintf('\n \t reflection')
                
                for z=1:domain.nozones %check in which zone is the point, to be checked
                    if inpolygon(yout(end-1,1), yout(end-1,2),  domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  )
                        if z==1
                            mesh=mesh1;
                        elseif z==2
                            mesh=mesh2;
                        end
                    end
                end
                
                [k1] = dsearchn(mesh,[y(end,1) y(end,2)]);
                [ktemp] = dsearchn([mesh(k1+1,:); mesh(k1-1,:)],[y(end,1) y(end,2)]);
                
                if ktemp==1
                    k2=k1+1;
                else
                    k2=k1-1;
                end
                
                V1=mesh(k1,:);
                V2=mesh(k2,:);
                
                P1=[yout(end-1,1) yout(end-1,2)];
                P2=[yout(end,1) yout(end,2)];
                u=P2-P1;
                v=V2-V1;
                
                alpha=atan(u(2)/u(1))*180/pi;
                if u(1)<0
                    alpha=180+alpha;
                end
                
                % there's difference in the eikonal angle
                % and the vector angle. why??
                % might be used directly alpha_chi
                
                %                 alpha_chi=atan(yout(end,4)/yout(end,3))*180/pi;
                %                 if yout(end,3)<0
                %                     alpha_chi=180+alpha_chi;
                %                 end
                %
                alpha_chi=atan(yout(end-1,4)/yout(end-1,3))*180/pi;
                if yout(end,3)<0
                    alpha_chi=180+alpha_chi;
                end
                
                % 2 different ways of computing reflection, can be tested
                beta=atan(v(2)/v(1))*180/pi; %this needs to be tested better
                alpha_reflection1=360-alpha+2*beta;
                
                
                if (alpha_reflection1>360)
                    alpha_reflection1=alpha_reflection1-360;
                end
                
                Ri=null(P2-P1);
                N=null(V2-V1);
                Rr = u - 2 * null(V2-V1)'* (dot(u,null(V2-V1)));
                
                
                alpha_reflection=atan(Rr(2)/Rr(1))*180/pi;
                
                
                
                if Rr(1)<0
                    alpha_reflection=180+alpha_reflection;
                end
                
                if alpha_reflection<0
                    alpha_reflection=360+alpha_reflection;
                end
                
                
                % if abs(alpha-alpha_chi)>0.1
                
                if abs(alpha-alpha_chi)>1
                    alpha_chi2=atan(yout(end,4)/yout(end,3))*180/pi;
                    if yout(end,3)<0
                        alpha_chi2=180+alpha_chi2;
                    end
                    fprintf('\n \t vector angle=%f,  eikonal angle (t-1)=%f, eikonal angle (t)=%f, ',alpha,alpha_chi,alpha_chi2);
                    
                end
                
                
                if abs(alpha_reflection-alpha_reflection1)>0.1
                    
                    fprintf('\n \t reflection angle=%f,  angle2 =%f , initial angle =%f',alpha_reflection,alpha_reflection1,itdir(1,a));
                    
                    break
                end
                
                %here could happen that cut-off is not recognized
                n_0 = refractive_interpolation(y(end,1),y(end,2),domain);
                
                y0(3)=n_0*cos(alpha_reflection*pi/180);
                y0(4)=n_0*sin(alpha_reflection*pi/180);
                
                y0(1) = y(end,1)+0.05*y0(3);
                y0(2) = y(end,2)+0.05*y0(4);
                
                % same values as last step
                y0(5) = y(end,5);
                y0(6) = y(end,6);
                y0(7) = y(end,7);
                
                if nt>refine
                    options = odeset(options,'InitialStep',tout(nt)-tout(nt-refine));
                end
                
                tstart = tout(end);
                
            elseif n_end<cut_off_limit && checkifinsidedomain([yout(end,1) yout(end,2)], domain)
                
                fprintf('\n \t cut off ');
                
                P1=[yout(end-1,1) yout(end-1,2)];
                P2=[yout(end,1) yout(end,2)];
                u=P2-P1;
                
                N=[nx_end ny_end]/(sqrt(nx_end^2+ny_end^2));
                
                Rr = u - 2 * N* (dot(u,N));
                
                
                alpha_reflection=atan(Rr(2)/Rr(1))*180/pi;
                
                if Rr(1)<0
                    alpha_reflection=180+alpha_reflection;
                end
                
                if alpha_reflection<0
                    alpha_reflection=360+alpha_reflection;
                end
                
                
                y0(3)=n_end*cos(alpha_reflection*pi/180);
                y0(4)=n_end*sin(alpha_reflection*pi/180);
                
                y0_3=cos(alpha_reflection*pi/180);
                y0_4=sin(alpha_reflection*pi/180);
                
                %to move away from cutoff conditions
                y0(1) = yout(end,1)+0.005*y0(3);
                y0(2) = yout(end,2)+0.005*y0(4);
                y0(5) = yout(end,5);
                y0(6) = yout(end,6);
                y0(7) = yout(end,7);
                
                
                options = odeset(options,'InitialStep',tout(nt)-tout(nt-refine));
                
                tstart = tout(nt);
                
            else
                
                flag_out=0;
                fprintf('\n \t outside')
                
            end
            
            
        end % while flagout
        
        itpo(1,1:length(yout(:,1)),a)= yout(:,1); % x
        itpo(2,1:length(yout(:,1)),a)= yout(:,2); % y
        
        itpo(3,1:length(yout(:,1)),a)= yout(:,7); %length of path
        
        % this can be done in different ways, can be ugraded.
        %1) with outputfcn in the ode
        %2) by simply calling refractive interpolation instead of eikonal
        % maybe these solutions are performing better, but for the moment works
        
        %!! this is very very slow
        %         for  i=1:length(yout(:,1))
        %             [~,itpo(4,i,a)]= eikonal_equation(tout(i).', yout(i,:).',domain); % refractive along path (not integrated)
        %         end
        
        % this should work properly
        for  i=1:length(yout(:,1))
            itpo(4,i,a)= sqrt(yout(i,3)^2+yout(i,4)^2); % refractive along path (not integrated)
        end
        
        % to be checked if we need k0 or it's already there
        itpo(5,1:length(yout(:,1)),a)= yout(:,6); %  absorption integrated
        itpo(6,1:length(yout(:,1)),a)= yout(:,5); % optical path length
        
        % to be fixed
        for  i=2:length(yout(:,1))
            P1=[yout(i-1,1) yout(i-1,2)];
            P2=[yout(i,1) yout(i,2)];
            u=P2-P1;
            itdir(i,a)=atan(u(2)/u(1))*180/pi;
            if u(1)<0
                itdir(i,a)=180+itdir(i,a);
            end
        end
        
      
    end % for max angles
    
end %% function raytracting



%pause


function [solution,n] = eikonal_equation(t,y,domain)
    
    
    xi=y(1);
    yi=y(2);
    zi=y(3);
    wi=y(4);
    
    
    
    [n, nx, ny, alpha] = refractive_interpolation(xi, yi, domain);
        
    yp = zeros(7,1);
    
    yp(1) = zi / (n^2) ;
    yp(2) = wi / (n^2) ;
    yp(3) = nx/n ;
    yp(4) = ny/n ;
    yp(5) = n;
    yp(6) = alpha;
    yp(7) = sqrt(yp(1)^2+yp(2)^2);
    
    solution = [yp(1); yp(2); yp(3); yp(4);  yp(5); yp(6); yp(7)];
    
    
    
end


function [n,nx,ny,alpha] = refractive_interpolation(x,y,domain)
    
    isInside=checkifinsidedomain([x y],domain);
    
    if isInside
        
        n=interpolation(domain,[x y],domain.nova-1);
        nx=interpolation(domain,[x y],domain.nova+1);
        ny=interpolation(domain,[x y],domain.nova+2);
        alpha=interpolation(domain,[x y],domain.nova);
        
    else
        
        %closest point
        k=0;
        k1=zeros(domain.nozones);
        dist=zeros(domain.nozones);
        
        for z=1:domain.nozones
            [k1(z), dist(z)]=dsearchn([ domain.( strcat('zone',num2str(z)) ).variables(1,:)' ...
                domain.( strcat('zone',num2str(z)) ).variables(2,:)' ],[x y]);
                        
        end
        
        if domain.nozones>1
            if dist(1)<dist(2)
                k=k1(1);
                z=1;
            else
                k=k1(2);
                z=2;
            end
        else
            k=k1(1);
            z=1;
        end
        x_close=domain.( strcat('zone',num2str(z)) ).variables(1,k);
        y_close=domain.( strcat('zone',num2str(z)) ).variables(2,k);
        
        n=interpolation(domain,[x_close y_close],domain.nova-1);
        nx=interpolation(domain,[x_close y_close],domain.nova+1);
        ny=interpolation(domain,[x_close y_close],domain.nova+2);
        alpha=interpolation(domain,[x_close y_close],domain.nova);
        
        
    end
    
end

function [value,isterminal,direction] = cutoff(t,y,domain)
    
    flag_out=checkifinsidedomain([y(1) y(2)],domain);
    
    if flag_out~=1
        
        value = flag_out;
        isterminal = 1;            % Stop at local minimum
        direction = 0;
        
    else
        
        n = refractive_interpolation(y(1), y(2), domain);
        flag=1;
        
        if n<0.001
            flag=0;
        end
        
        value = flag;
        isterminal = 1;            % Stop at local minimum
        direction = 0;            % [local minimum, local maximum]
        
    end
end
