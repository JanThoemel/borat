function value=interpolation(domain,po,index)
%% this function returns an interpolated value

%%  input
%%
%%  output
%%
%%  internal parameters
    value=0;     

    %% find triangle in which the po is inside
    zonenumber=0;
    for z=1:domain.nozones
        %!checkuseofvariables!
        if inpolygon(po(1),po(2),  domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  )  
            zonenumber=z;
        end
    end
    
    if zonenumber
        %% ti = pointLocation(T,QP) returns the enclosing triangles or tetrahedra for each query point in QP. Each row in matrix QP contains the x, y, (and possibly z) coordinates of a query point.
        triangleindex=pointLocation(     domain.( strcat('zone',num2str(zonenumber)) ).DT  ,  po(1),po(2)   );
        if isfinite(triangleindex)
        
        pointindices=domain.( strcat('zone',num2str(zonenumber)) ).DT.ConnectivityList( triangleindex,: );

        triangle(:,1)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(1));
        triangle(:,2)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(2));
        triangle(:,3)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(3));

        value= griddata(triangle(1,:),triangle(2,:),triangle(index,:),po(1),po(2),'cubic');  
        else
            return
    end
end %% function interpolation
