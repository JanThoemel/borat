function zonenumber=checkifinsidedomain3D(po,domain)
    %%  internal parameters
    zonenumber=0;
    
    %% find closest triangle i.e. its coordinates and indices
    %% convert mesh to vector
    for z=1:domain.nozones
        
        %             if inpolyhedron(domain.( strcat('zone',num2str(z)) ).BoundaryFaces, domain.( strcat('zone',num2str(z)) ).BoundaryPoints, po(1),po(2), po(3))
        %                 zonenumber=1;
        %             end
        
        if isinf(po(1))
            po(1)
        end
        
        if isinf(po(2))
            po(2)
        end
        if isinf(po(3))
            po(3)
        end
        
        tetraindex=pointLocation(domain.( strcat('zone',num2str(z)) ).DT, po(1),po(2), po(3));
        
        if isfinite(tetraindex)
            
            zonenumber=1;
            
        end
        
    end
    
end
