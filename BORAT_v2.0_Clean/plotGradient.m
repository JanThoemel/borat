figure
    hold on
    
    for z=1:domain.nozones
        
        [domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,:),domain.( strcat('zone',num2str(z)) ).variables(domain.nova+2,:)] = trigradient(domain.( strcat('zone',num2str(z)) ).variables(1,:)',domain.( strcat('zone',num2str(z)) ).variables(2,:)',domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)',domain.( strcat('zone',num2str(z)) ).delaunay);
        
    end
    
    for z=1:domain.nozones
        %%plot boundary
        plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)
        %% plot coll freq
        %h=trisurf(domain.(strcat('zone',num2str(z))).delaunay,squeeze(domain.(strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)));
        %  [CS,h]=tricontf(squeeze(domain.(strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,domain.(strcat('zone',num2str(z))).delaunay,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)));
        %[CS,h]=TRICONTF(X,Y,M,Z)
        %h=trisurf(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)));
        %[C,h]=tricontour(domain.(strcat('zone',num2str(z))).delaunay,squeeze(domain.(strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)),nb_lines);
        %VC=[0,0.0001,0.001,0.01,0.10,1];
        
        grad=sqrt(domain.( strcat('zone',num2str(z)) ).variables(domain.nova+1,:).^2+ domain.( strcat('zone',num2str(z)) ).variables(domain.nova+2,:).^2);
        
        
        h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(grad) );
        set(h, 'edgecolor','none');
        axis tight;caxis([0 50]);

        hcb=colorbar;title(hcb,'gradient                 ');
        box on;grid off;
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
        
    end
    hold off