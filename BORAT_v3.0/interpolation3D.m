function value=interpolation3D(domain,po,index)
    % this function returns an interpolated value
    
    value=0;
    
    % find triangle in which the po is inside
    zonenumber=0;
    for z=1:domain.nozones
        
        
        tetraindex=pointLocation(domain.( strcat('zone',num2str(z)) ).DT, po(1),po(2), po(3));
        
        if isfinite(tetraindex)
            
            zonenumber=z;
            
        end
    end
    
    if zonenumber
        % ti = pointLocation(T,QP) returns the enclosing triangles or tetrahedra for each query point in QP. Each row in matrix QP contains the x, y, (and possibly z) coordinates of a query point.
%         tetraindex=pointLocation(     domain.( strcat('zone',num2str(zonenumber)) ).DT  ,  po(1),po(2), po(3)   );
        pointindices=domain.( strcat('zone',num2str(zonenumber)) ).DT.ConnectivityList( tetraindex,: );
        
        tetra(:,1)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(1));
        tetra(:,2)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(2));
        tetra(:,3)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(3));
        tetra(:,4)=domain.( strcat('zone',num2str(zonenumber)) ).variables(:,pointindices(4));
        
        value= griddata(tetra(1,:),tetra(2,:),tetra(3,:),tetra(index,:),po(1),po(2),po(3),'natural');
    end
end %% function interpolation


