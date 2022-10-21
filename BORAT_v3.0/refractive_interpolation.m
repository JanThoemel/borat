function [n,nx,ny,nz] = refractive_interpolation(x,y,z,domain)
    
    isInside=checkifinsidedomain3D([x y z],domain);
    
    if isInside
        
        n=interpolation3D(domain,[x y z],domain.nova-3);
        nx=interpolation3D(domain,[x y z],domain.nova-2);
        ny=interpolation3D(domain,[x y z],domain.nova-1);
        nz=interpolation3D(domain,[x y z],domain.nova);
        
        if n<0.001
            n=0.00;
        end
    else
        
        n=1;
        nx=0;
        ny=0;
        nz=0;
        
    end
        
    
end