function solution = eikonal3D(t,y,domain)
    
    
    xi=y(1);
    yi=y(2);
    zi=y(3);
    ui=y(4);
    vi=y(5);
    wi=y(6);
    
    if isinf(xi)
        xi;
    end
    
    if isnan(xi)
        xi;
    end
    
    [n, nx, ny,nz] = refractive_interpolation(xi, yi, zi,domain);
    
    if n==0
        n;
        n=0.009;
    end
    
    if n<0.01
        n;
    end
    
    if isinf(n)
        n;
    end
    
    if isnan(n)
        n;
    end
    
    yp = zeros(8,1);
    
    yp(1) = ui / (n^2) ;
    yp(2) = vi / (n^2) ;
    yp(3) = wi / (n^2) ;
    
    yp(4) = nx/n ;
    yp(5) = ny/n ;
    yp(6) = nz/n ;
    yp(7) = n;
    yp(8) = sqrt(yp(1)^2+yp(2)^2+yp(3)^2); %path_length
    
    solution = [yp(1); yp(2); yp(3); yp(4); yp(5); yp(6); yp(7); yp(8)];
    
end