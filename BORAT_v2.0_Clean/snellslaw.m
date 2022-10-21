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
