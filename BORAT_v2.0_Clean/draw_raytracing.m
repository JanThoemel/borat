function draw_raytracing(domain,plotscdir,scdir,scdir_range,collisionmodel,absorptionlimits,pooo1,pooo2,itpo,maxangles,itdir,symmetrylineencounter,writefig,runFolder)
    
    cd (runFolder);
    
    
    % compute direction and range to space craft
    if plotscdir
        if ( pooo1(1)==pooo2(1) && pooo1(2)==pooo2(2) ) %% compute and show spacecraft direction and range only if raytracing start from antenna location i.e. pooo1=pooo2
            for it=1:3
                scdir_x(it)=pooo1(1)+cos(scdir/180*pi)*(it-1)*3;                 scdir_y(it)=pooo1(2)+sin(scdir/180*pi)*(it-1)*3;
                scdirr1_x(it)=pooo1(1)+cos((scdir+scdir_range)/180*pi)*(it-1)*3; scdirr1_y(it)=pooo1(2)+sin((scdir+scdir_range)/180*pi)*(it-1)*3;
                scdirr2_x(it)=pooo1(1)+cos((scdir-scdir_range)/180*pi)*(it-1)*3;	scdirr2_y(it)=pooo1(2)+sin((scdir-scdir_range)/180*pi)*(it-1)*3;
            end
        end
    end
    % plot rays and spacecraft direction/direction range within refractive index contour plot
    figure
    hold on;
    V=[1e-17,1e-16,1e-15,1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-3:0.02:1-1e-3,1-1e-3,1-1e-4,1-1e-5,1-1e-6,1-1e-7,1-1e-8,1-1e-9,1-1e-10,1-1e-11,1-1e-12,1-1e-13,1-1e-14,1-1e-15,1-1e-16,1-1e-17,1];
    %V=[0:0.1:0.9,0.9:0.01:1];
    %V=[0.09:0.01:0.99];
    %axis equal;
    for z=1:domain.nozones
        %%plot flowfield with ri, raytraces and s/c direction/range
        
        %%plot ri contour plot
        if 0
            %% plot boundary
            plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)
            %% plot ri
            [C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),V);
        end
        %%plot ri surface plot
        if 1
            %% plot boundary
            plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)
            %% plot ri
            h=trisurf(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)));
            set(h, 'edgecolor','none');
        end
        maxX(z) = max(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        minX(z) = min(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        maxY(z) = max(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        minY(z) = min(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        maxC(z) = max(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)  );
        minC(z) = min(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)  );
        maxX = max(maxX);
        minX = min(minX);
        maxY = max(maxY);
        minY = min(minY);
        maxC = max(maxC);
        minC = min(minC);
    end
    
    z_max = max(max(get(h,'Zdata')));
    scdir_z=[z_max,z_max,z_max];
    
    %%plot rays in domain
    itpoz=z_max*ones(   size( itpo(1,:,1),2)   ,1);%z coordinates of lines to be plotted
    for a=1:maxangles
        %!can they be colored according to attenuation?
        plot3(  itpo(1,1:nnz(itpo(1,:,a)),a)  ,  itpo(2,1:nnz(itpo(1,:,a)),a),itpoz(1:nnz(itpo(1,:,a))),'k-','LineWidth',1,'color',[105/256 105/256 105/256]);
    end
    %%plot s/c direction and direction range
    if plotscdir
        if ( pooo1(1)==pooo2(1) && pooo1(2)==pooo2(2) )%% compute and show spacecraft direction and range only if raytracing start from antenna location i.e. pooo1=pooo2
            plot3(scdir_x,scdir_y,scdir_z,'k-.','color',[20/256 20/256 20/256]);
            %plot3(scdirr1_x,scdirr1_y,scdir_z,'k--','color',[20/256 20/256 20/256]);
            %plot3(scdirr2_x,scdirr2_y,scdir_z,'k--','color',[20/256 20/256 20/256]);
        end
    end
    text(1.5,3.1,z_max+1,'direction to TGO');
    set(h, 'edgecolor','none');
    set(gcf, 'Position', [10, 10, 1100, 800]);
    % xlim([-0.17 6]);ylim([0 5.185]);caxis([0 1])
        
         axis equal; xlim([minX maxX]); ylim ([minY maxY]); caxis([minC maxC]);
    
    hcb=colorbar;title(hcb,'\mu');
    box off;grid off;
    set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
    hold off;
    if writefig
        savefig('raytracing')
        %print('-painters','raytracing','-depsc')
        print('-painters','raytracing','-dpng')
    end
    
    
    figure
    hold on
    for a=1:maxangles
        plot(itpo(3,1:nnz(itpo(1,:,a)),a),itpo(4,1:nnz(itpo(1,:,a)),a))
    end
    xlabel('path length [m]');ylabel('refractive index [-]');axis([0 inf -0.05 1.05]);
    hold off
    if writefig
        savefig('path_mu')
        %print('-painters','path_mu','-depsc')
        print('-painters','path_mu','-dpng')
    end
    %% plot attenuation over path length
    if collisionmodel~=0 %% no plot if collisions/absorption was not accounted for
        figure
        eps = 2e-16;
        for a=1:maxangles
            semilogy(itpo(3,1:nnz(itpo(1,:,a)),a), max(eps, itpo(5,1:nnz(itpo(1,:,a)),a)));
            hold on
        end
        xlabel('path length [m]');ylabel('attenutation [dB]');axis([0 inf absorptionlimits]);
        %% absorption cut off
        semilogy([0 inf], [absorptionlimits(2)-10 absorptionlimits(2)-10]);
        if writefig
            savefig('path_kappa')
            %print('-painters','path_kappa','-depsc')
            print('-painters','path_kappa','-dpng')
        end
        hold off
    end
    
    cd ..;
    
end