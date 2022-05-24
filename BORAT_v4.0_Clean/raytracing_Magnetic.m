function [itdir, itpo,symmetrylineencounter,Y_int_1,Y_int_2,Y_int_3,YY]=raytracing_Magnetic(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, absorptioncutoff,symmetryline,ss,f,Cartesian,B1,B2,B3,MagneticField)
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
    for z=1:domain.nozones
        
        [domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,:),domain.( strcat('zone',num2str(z)) ).variables(domain.nova+2,:)] = trigradient(domain.( strcat('zone',num2str(z)) ).variables(1,:)',domain.( strcat('zone',num2str(z)) ).variables(2,:)',domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)',domain.( strcat('zone',num2str(z)) ).delaunay);
        %% parameters to calculate magnetic refractive index
        %% speed of light, wikipedia 8/7/2017, m/s
         c=299792458;
        %% kk=e^2/(4 pi^2 e0 m)
        kk=80.5;
        %% calculated in readflowfield_tecplot.m
        electronnumberdensity = domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:);
        %% ee = e/2*pi*m_e
        ee= (1.602176634e-19)/(2*pi*f*(9.1093837015e-31));
        %% X parameter of ri
        XX=kk*electronnumberdensity/f^2;
    end

    %% Get B-field from File --> additional work to be done for external B-field
    
    if Cartesian == 0 %Array of B-field vectors in cartesian coords
        for i=1:domain.(strcat('zone',num2str(z))).N
            [B_vectors(1,i),B_vectors(2,i),B_vectors(3,i)]=pol2cart(domain.(strcat('zone',num2str(z))).variables(B2,i),domain.(strcat('zone',num2str(z))).variables(B1,i),domain.(strcat('zone',num2str(z))).variables(B3,i));
        end
    elseif Cartesian == 1
        B_vectors=[domain.(strcat('zone',num2str(z))).variables(B1,:);domain.(strcat('zone',num2str(z))).variables(B2,:);domain.(strcat('zone',num2str(z))).variables(B3,:)];
    end
    %B-field quiver for vector field
    B_vectors(3,:)=0;

    domain.(strcat('zone',num2str(z))).variables(domain.nova+9,:)=B_vectors(1,:);
    domain.(strcat('zone',num2str(z))).variables(domain.nova+10,:)=B_vectors(2,:);
    domain.(strcat('zone',num2str(z))).variables(domain.nova+11,:)=B_vectors(3,:);
    
    %magnitude of B-field to calculate ri
    for i=1:length(B_vectors)
        B_vec_mag(i)= sqrt(B_vectors(1,i)^2 +B_vectors(2,i)^2 +B_vectors(3,i)^2);
    end
    Y=ee*B_vec_mag;
    %% iteration points: 1: x, 2: y, 3: pathlength, 4: ri, 5: ec
    itpo=zeros(5,maxsteps,maxangles);
    %itpo_Magnetic=zeros(5,maxsteps,maxangles);
    %% iteration directions
    itdir=zeros(maxsteps,maxangles);
    %itdir_Magnetic=zeros(maxsteps,maxangles);

    %% starting directions and points for all rays
    if maxangles==1
        itdir(1,1)=dir1;
        itpo(1,1,1)=pooo1(1);
        itpo(2,1,1)=pooo1(2);


        %calculate ray vector from itpo and itdir
        %|k_vec|=2*pi*f/c
        RayVector= [(2*pi*f/c)*cos(itdir(1,1)),(2*pi*f/c)*sin(itdir(1,1)),0]; %%only 2D rays
        
        %Find B-field vector at starting point
        Y_int_1(1,1,1) = interpolation(domain,itpo(1:2,1,1)',domain.nova+9);
        Y_int_2(1,1,1) = interpolation(domain,itpo(1:2,1,1)',domain.nova+10);
        Y_int_3(1,1,1) = interpolation(domain,itpo(1:2,1,1)',domain.nova+11);
        B_Vector_at_itpo = [Y_int_1,Y_int_2,Y_int_3];
        if norm(B_Vector_at_itpo) == 0
                    Y(1,1)=0;
                    YY(1,1,1)=0;
                    Theta(1,1)=0;
                    
                    if 1-XX>=0
                        %%refractive index
                        ri=sqrt(1-XX);
                        %%absorption coefficient
                        ac=0;
                    elseif 1-XX<0 %is this condition physical?
                        warning('ri_ac_domain_100')
                        %%refractive index
                        ri=0;
                        %%absorption coefficient, check sign convention
                        ac=2*pi*f/c*sqrt(XX-1);
                        domain.( strcat('zone',num2str(z)) ).variables(domain.nova+12,:)=ac;
                        domain.( strcat('zone',num2str(z)) ).variables(domain.nova+13,:)=ac;
                    end
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+3,:)=ri;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+4,:)=ri;
          else
                   
                    Theta=atan2(norm(cross(RayVector,B_Vector_at_itpo)), dot(RayVector,B_Vector_at_itpo));
        
                    %refractive index at start point
                    
                    for i=1:length(electronnumberdensity)
                            Y(i)=ee*B_vec_mag(i);
                           % Y(i)=0;
                            ri_1(i)=((1-XX(i))^2)*((Y(i))^2)*(cos(Theta(1,1))^2);
                            ri_2(i)=(0.5*((Y(i))^2)*(sin(Theta(1,1))^2))^2;
                            ri_3(i)=sqrt(ri_2(i)+ri_1(i));
                            ri_4(i)=1-XX(i)-(0.5*((Y(i))^2)*(sin(Theta(1,1))^2));

                            ri_5(i)=ri_4(i)+ri_3(i);
                            ri_6(i)=ri_4(i)-ri_3(i);

                            ri_7(i)=XX(i)*(1-XX(i));

                            ri_8(i)=1-(ri_7(i)/ri_5(i));
                            ri_9(i)=1-(ri_7(i)/ri_6(i));

                            riPlus(i)=sqrt(ri_8(i));
                            riMinus(i)=sqrt(ri_9(i));
                            ac_Plus(i)=0;
                            ac_Minus(i)=0;
                            if isreal(riPlus(i))==false
                                riPlus(i)=0;
                                ac_Plus(i)=2*pi*f/c*sqrt(XX(i)-1);
                            elseif isreal(riMinus(i))==false
                                riMinus(i)=0;
                                ac_Minus(i)=2*pi*f/c*sqrt(XX(i)-1);
                            end
                    
                    end
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+14,:)=Y;
                    YY(1,1,1)=interpolation(domain,itpo(:,1,1),domain.nova+14);
                  
                    riPlus(isinf(riPlus) | isnan(riPlus))=1;
                    riMinus(isinf(riMinus) | isnan(riMinus))=1;  
                    
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+3,:)=riPlus;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+4,:)=riMinus;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+12,:)=ac_Plus;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+13,:)=ac_Minus;
            end
        
        %% starting point refractive index for all rays
        itpo(4,1,:)=interpolation(domain,itpo(:,1,1),domain.nova+3);
        %% starting point absorption coefficient for all rays
        itpo(5,1,:)=interpolation(domain,itpo(:,1,1),domain.nova+12);
    else
        for a=1:maxangles
            %% starting direction (fan)
            itdir(1,a)=dir1-anglestep*(a-1);
            %% initialpoints
            itpo(1,1,a)=pooo1(1)-x_step*(a-1);
            itpo(2,1,a)=pooo1(2)-y_step*(a-1);
            %calculate ray vector from itpo and itdir
            %|k_vec|=2*pi*f/c
            RayVector= [(2*pi*f/c)*cos(itdir(1,a)),(2*pi*f/c)*sin(itdir(1,a)),0]; %%only 2D rays
            
            %Find B-field vector at starting point
            Y_int_1(1,1,a) = interpolation(domain,itpo(1:2,1,a)',domain.nova+9);
            Y_int_2(1,1,a) = interpolation(domain,itpo(1:2,1,a)',domain.nova+10);
            Y_int_3(1,1,a) = interpolation(domain,itpo(1:2,1,a)',domain.nova+11);
            B_Vector_at_itpo = [Y_int_1(1,1,a),Y_int_2(1,1,a),Y_int_3(1,1,a)];

            if norm(B_Vector_at_itpo) == 0
                    
                    YY(1,1,a)=0;
                    Theta(1,a)=0;
                    
                    if 1-XX>=0
                        %%refractive index
                        ri=sqrt(1-XX);
                        %%absorption coefficient
                        ac=0;
                    elseif 1-XX<0 %is this condition physical?
                        warning('ri_ac_domain_100')
                        %%refractive index
                        ri=0;
                        %%absorption coefficient, check sign convention
                        ac=2*pi*f/c*sqrt(XX-1);
                    end
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+3,:)=ri;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+4,:)=ri;
            else
                    Theta(1,a)=atan2(norm(cross(RayVector,B_Vector_at_itpo)), dot(RayVector,B_Vector_at_itpo));
                    
                    
                    %refractive index at start point
                    
                    for i=1:length(electronnumberdensity)
                            Y(i)=ee*B_vec_mag(i);
                         
                            ri_1(i)=((1-XX(i))^2)*((Y(i))^2)*(cos(Theta(1,a))^2);
                            ri_2(i)=(0.5*((Y(i))^2)*(sin(Theta(1,a))^2))^2;
                            ri_3(i)=sqrt(ri_2(i)+ri_1(i));
                            ri_4(i)=1-XX(i)-(0.5*((Y(i))^2)*(sin(Theta(1,a))^2));

                            ri_5(i)=ri_4(i)+ri_3(i);
                            ri_6(i)=ri_4(i)-ri_3(i);

                            ri_7(i)=XX(i)*(1-XX(i));

                            ri_8(i)=1-(ri_7(i)/ri_5(i));
                            ri_9(i)=1-(ri_7(i)/ri_6(i));

                            riPlus(i)=sqrt(ri_8(i));
                            riMinus(i)=sqrt(ri_9(i));
                            ac_Plus(i)=0;
                            ac_Minus(i)=0;
                            if isreal(riPlus(i))==false
                                riPlus(i)=0;%real(riPlus(i));
                                ac_Plus(i)=2*pi*f/c*sqrt(XX(i)-1);
                            elseif isreal(riMinus(i))==false
                                riMinus(i)=0;
                                ac_Minus(i)=2*pi*f/c*sqrt(XX(i)-1);
                            end
                    
                    end
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+14,:)=Y;
                    YY(1,1,a)=interpolation(domain,itpo(:,1,1),domain.nova+14);
                  
                    riPlus(isinf(riPlus) | isnan(riPlus))=1;
                    riMinus(isinf(riMinus) | isnan(riMinus))=1;  
                    
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+3,:)=riPlus;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+4,:)=riMinus;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+12,:)=ac_Plus;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+13,:)=ac_Minus;
                end
            %% starting point refractive index for all rays
            itpo(4,1,a)=interpolation(domain,itpo(:,1,a),domain.nova+3);
            %% starting point absorption coefficient for all rays
            itpo(5,1,a)=interpolation(domain,itpo(:,1,a),domain.nova+12);
        end %% for
    end %%if    

    fprintf('\n  starting point refractive index for all rays %.6f',  itpo(4,1,1));
    fprintf('\n  starting point absorption coefficient for all rays %.6f ',itpo(5,1,1))
    
    %mesh1 and mesh2 contains the point in the wall surface for each
    %zone
    %for now it works in 2D with 2 zones symmetrical
    
    
    mesh1=[domain.zone1.variables(1,domain.zone1.bound(:))',domain.zone1.variables(2,domain.zone1.bound(:))'];
    
    if (domain.nozones==2)
        mesh2=[domain.zone2.variables(1,domain.zone2.bound(:))',domain.zone2.variables(2,domain.zone2.bound(:))'];
    end
    


    for a=1:maxangles
        s=1;msg=' ';
        %% while loop to check if ray-tracing stepping leaves 1) domain 2) maxsteps are exceeded or attenuation is above threshold
        %% if so go to next ray
        %! evaluate and improve efficiency of this loop
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while checkifinsidedomain(itpo(1:2,s,a)',domain)==1 && s<maxsteps  && itpo(5,s,a)<absorptioncutoff
            
            boundary=0;
            s=s+1;
            if s > 2
                %calculate ray vector from itpo and itdir
                %|k_vec|=2*pi*f/c
                RayVector= [(2*pi*f/c)*cos(itdir(s-1,a)),(2*pi*f/c)*sin(itdir(s-1,a)),0]; %%only 2D rays
                
                %Find B-field vector at starting point
                Y_int_1(1,s-1,a) = interpolation(domain,itpo(:,s-1,a)',domain.nova+9);
                Y_int_2(1,s-1,a) = interpolation(domain,itpo(:,s-1,a)',domain.nova+10);
                Y_int_3(1,s-1,a) = interpolation(domain,itpo(:,s-1,a)',domain.nova+11);
  
                B_Vector_at_itpo = [Y_int_1(1,s-1,a),Y_int_2(1,s-1,a),Y_int_3(1,s-1,a)];
                if norm(B_Vector_at_itpo) == 0
                    Y(s-1,a)=0;
                    Theta(s-1,a)=0;
                    
                    if 1-XX>=0
                        %%refractive index
                        ri=sqrt(1-XX);
                        %%absorption coefficient
                        ac=0;
                    elseif 1-XX<0 %is this condition physical?
                        warning('ri_ac_domain_100')
                        %%refractive index
                        ri=0;
                        %%absorption coefficient, check sign convention
                        ac=2*pi*f/c*sqrt(XX-1);
                    end
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+3,:)=ri;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+4,:)=ri;
                else
                    
                    Theta(s-1,a)=atan2(norm(cross(RayVector,B_Vector_at_itpo)), dot(RayVector,B_Vector_at_itpo));
                    
                    
                    %refractive index at start point

                    for i=1:length(electronnumberdensity)
                            Y(i)=ee*B_vec_mag(i);
                           
                            ri_1(i)=((1-XX(i))^2)*((Y(i))^2)*(cos(Theta(s-1,a))^2);
                            ri_2(i)=(0.5*((Y(i))^2)*(sin(Theta(s-1,a))^2))^2;
                            ri_3(i)=sqrt(ri_2(i)+ri_1(i));
                            ri_4(i)=1-XX(i)-(0.5*((Y(i))^2)*(sin(Theta(s-1,a))^2));

                            ri_5(i)=ri_4(i)+ri_3(i);
                            ri_6(i)=ri_4(i)-ri_3(i);

                            ri_7(i)=XX(i)*(1-XX(i));

                            ri_8(i)=1-(ri_7(i)/ri_5(i));
                            ri_9(i)=1-(ri_7(i)/ri_6(i));

                            riPlus(i)=sqrt(ri_8(i));
                            riMinus(i)=sqrt(ri_9(i));
                            ac_Plus(i)=0;
                            ac_Minus(i)=0;
                            if isreal(riPlus(i))==false
                                riPlus(i)=0;%real(riPlus(i));
                                ac_Plus(i)=2*pi*f/c*sqrt(XX(i)-1);
                            elseif isreal(riMinus(i))==false
                                riMinus=0;
                                ac_Minus(i)=2*pi*f/c*sqrt(XX(i)-1);
                            end
        
                    end
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+14,:)=Y;
                    YY(1,s-1,a)=interpolation(domain,itpo(:,1,1),domain.nova+14);
                  
                    riPlus(isinf(riPlus) | isnan(riPlus))=1;
                    riMinus(isinf(riMinus) | isnan(riMinus))=1;  
                    
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+3,:)=riPlus;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+4,:)=riMinus;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+12,:)=ac_Plus;
                    domain.( strcat('zone',num2str(z)) ).variables(domain.nova+13,:)=ac_Minus;
                end
            end

            po_far(1)=itpo(1,s-1,a)+cos(itdir(s-1,a)/180*pi)*1.5*ss;
            po_far(2)=itpo(2,s-1,a)+sin(itdir(s-1,a)/180*pi)*1.5*ss;
            %                    if not(checkifinsidedomain(po_far,domain))
            %! check this logic!
            if not(  checkifinsidedomain(po_far,domain)  ) && not(po_far(2) < 0 && symmetryline)
                
                %check if po_far hits the body surface
                if not (checkifboundary(po_far,domain))
                    
                    fprintf('\tpo_far outside')
                    break;
                else
                    fprintf('\t reflection')
                    boundary=1;
                end
            end
            %if po_far outside any zone - break raytracing of this ray
            
            if boundary %% improve algorithm of reflection
                
                P1=itpo(1:2,s-1,a)';
                
                X=[itpo(1,s-1,a) po_far(1)];
                Y=[itpo(2,s-1,a) po_far(2)];
                
                if(P1(2)>0)
                    mesh=mesh1;
                else
                    mesh=mesh2;
                end
                
                [xi,yi] = polyxpoly(X,Y,mesh(:,1),mesh(:,2));
                
                [k1] = dsearchn(mesh,[xi yi]);
                [ktemp] = dsearchn([mesh(k1+1,:); mesh(k1-1,:)],[xi yi]);
                
                if ktemp==1
                    k2=k1+1;
                else
                    k2=k1-1;
                end
                
                V1=mesh(k1,:);
                V2=mesh(k2,:);
                
                P_inter=[xi yi];
                u=P_inter-P1;
                v=V2-V1;
                beta=atan(v(2)/v(1))*180/pi; %this needs to be tested better
                
                
                if(itdir(s-1,a)>180)
                    alpha_reflection=360-itdir(s-1,a)+2*beta;
                    %                     else %do we need another option?
                    %                         alpha_reflection=360-itdir(s-1,a)+2*beta;
                end
                
                if (alpha_reflection>360)
                    alpha_reflection=alpha_reflection-360;
                end
                
                fprintf('\n \treflection angle=%f, initial angle =%f , body angle =%f',alpha_reflection,itdir(s-1,a),atan(v(2)/v(1))*180/pi);
                
                %%%%%%
                
                %once reflection is found, set s-step to intersection, and s+1 bounced
                itpo(1,s,a)=xi;
                itpo(2,s,a)=yi;
                itpo(3,s,a)=itpo(3,s-1,a)+sqrt(  (itpo(1,s,a)-itpo(1,s-1,a))^2 + (itpo(2,s,a)-itpo(2,s-1,a))^2   );
                itpo(4,s,a)=interpolation(domain,itpo(1:2,s,a)',domain.nova+3);
                itpo(5,s,a)= itpo(5,s-1,a)+interpolation(domain,itpo(1:2,s,a)',domain.nova+12)*(itpo(3,s,a)-itpo(3,s-1,a));
                itdir(s,a)=itdir(s-1,a);
                
                %assign next step to the bounce
                s=s+1;
                itpo(1,s,a)=xi+cos(alpha_reflection/180*pi)*ss/2;
                itpo(2,s,a)=yi+sin(alpha_reflection/180*pi)*ss/2;
                itpo(3,s,a)=itpo(3,s-1,a)+sqrt(  (itpo(1,s,a)-itpo(1,s-1,a))^2 + (itpo(2,s,a)-itpo(2,s-1,a))^2   );
                itpo(4,s,a)=interpolation(domain,itpo(1:2,s,a)',domain.nova+3);
                itpo(5,s,a)= itpo(5,s-1,a)+interpolation(domain,itpo(1:2,s,a)',domain.nova+12)*(itpo(3,s,a)-itpo(3,s-1,a));
                itdir(s,a)=alpha_reflection;
                
                continue %exit from this cycle
                
            end
            
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
            if (po_temp(2)<0 && symmetryline)
                itpo(1,s,a)=itpo(1,s-1,a);
                itpo(2,s,a)=itpo(2,s-1,a);
                itpo(3,s,a)=itpo(3,s-1,a)+sqrt(  (itpo(1,s,a)-itpo(1,s-1,a))^2 + (itpo(2,s,a)-itpo(2,s-1,a))^2   );
                itpo(4,s,a)=interpolation(domain,itpo(1:2,s,a)',domain.nova+3);
                %% integrate absorption coefficients
                itpo(5,s,a)= itpo(5,s-1,a)+interpolation(domain,itpo(1:2,s,a)',domain.nova+12)*(itpo(3,s,a)-itpo(3,s-1,a));
                itdir(s,a)=-itdir(s-1,a);
                symmetrylineencounter(a)=symmetrylineencounter(a)+1;
                continue
            end
            
            %!if po_far in next zone and po_temp in current zone, shorten stepsize
            %!if po_far and po_temp in next zone enlarge step size
            
            
            ri_po_temp=interpolation(domain,po_temp,domain.nova+3);
            
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
            itpo(4,s,a)=interpolation(domain,itpo(1:2,s,a)',domain.nova+3);
            %% integrate absorption coefficients over path
            itpo(5,s,a)=itpo(5,s-1,a)+interpolation(domain,itpo(1:2,s,a)',domain.nova+12)*(itpo(3,s,a)-itpo(3,s-1,a));
            
            
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
