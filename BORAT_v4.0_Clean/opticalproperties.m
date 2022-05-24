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
    %% Avogadro constant, wikipedia, mol-1
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
            Z=0;
            %%compute optical properties
            if 1-X>=0
                %%refractive index
                ri=sqrt(1-X);
                %%absorption coefficient
                ac=0;
            elseif 1-X<0 %is this condition physical?
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
                %summand1=summand1*0.999999;
                
                % fix of the sign 05/2021
                if summand1<0
                    summand1=-root1; %modified by vincent - 12/2020
                else
                    summand1=root1;
                end
                %pause
            end
            ri=         sqrt(0.5* (summand1+root1));
            ac=2*pi*f/c*sqrt(0.5* (-summand1+root1));
            
            if Z==0 %numerical bug that makes ac not equal to zero in some cases
                ac=0;
            end
            
            
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
