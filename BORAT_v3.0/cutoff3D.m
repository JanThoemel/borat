function [value,isterminal,direction] = cutoff3D(t,y,domain)
    
    flag_out=checkifinsidedomain3D([y(1) y(2) y(3)],domain);
    
    if flag_out~=1
        
        value = flag_out;
        isterminal = 1;            % Stop at local minimum
        direction = 0;
        
    else
        
        [n, nx, ny ,nz] = refractive_interpolation(y(1), y(2), y(3),domain);
        
        flag_cutoff=1;
        
        if n<0.001
            flag_cutoff=0;
        end
        
        value = flag_cutoff;
        isterminal = 1;            % Stop at local minimum
        direction = 0;            % [local minimum, local maximum]
        
    end
end