function boundary=checkifboundary(po,domain)

    %this simple function checks if the po(1,2) coordinates are 
    %in the range of the wall surface, delimited by limitx, limity
    
    %can be improved in the future, using domain to find the boundary
    
    boundary=0;

    if ( po(1)>domain.limitx1 && po(1)< domain.limitx2 && po(2)>domain.limity1 && po(2)<domain.limity2 )
        boundary=1;
    end
    
end
