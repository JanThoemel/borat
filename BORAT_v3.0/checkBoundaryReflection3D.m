function boundary=checkBoundaryReflection3D(po,domain)
    
    %inpolyhedron(domain.( strcat('zone',num2str(z)) ).Wall.BoundaryFaces, domain.( strcat('zone',num2str(z)) ).Wall.BoundaryPoints, po(1),po(2), po(3))
    
    % heavytest n>1 - to be tested
    
    boundary=0;
    for z=1:domain.nozones
        
        if intriangulation(domain.( strcat('zone',num2str(z)) ).Wall.variables(:,:)',domain.( strcat('zone',num2str(z)) ).Wall.NodeMap(:,:),po(:,:),50)
%             if ( po(1)>domain.limitx1 && po(1)< domain.limitx2 && po(2)>domain.limity1 && po(2)<domain.limity2  && po(3)>domain.limitz1 && po(3)<domain.limitz2)
                
                boundary=1;
                
%             end
        end
        
    end
    
end
