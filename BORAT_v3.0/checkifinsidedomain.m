function zonenumber=checkifinsidedomain(po,domain)
    %%  this function returns the zone and triangle (through the triangle's 3
    %%  vertex coordinates) in which the point is located. Returns zonenumber=0 if the point is outside the domain
    
    %%  input:
    %%  double(2)   po          coordinates of point that is checked to be inside domain (which zone and which triangle)
    %%  struct      domain      contains all flowfield variables
    %%  output:
    %%  int         zonenumber  number of zone in which the point po is located
    %%  double(2,3) triangle    coordinates of the three vertices of the triangle in which the point po is located
    %%  internal parameters
    zonenumber=0;
    
    %% find closest triangle i.e. its coordinates and indices
    %% convert mesh to vector
    for z=1:domain.nozones
        if inpolygon(po(1),po(2),  domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  )
            
            triangleindex=pointLocation(     domain.( strcat('zone',num2str(z)) ).DT  ,  po(1),po(2)   );
            if isfinite(triangleindex)
                %this has been changed due to errors on some cases, but
                %looks like it's really heavy in terms of computational
                %time
                zonenumber=1;
            end
        end
    end
