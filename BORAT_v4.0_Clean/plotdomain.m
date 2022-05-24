function plotdomain(domain,criticaldensity,max_edensity,disp_edensity,disp_coll_freq,disp_mu,disp_kappa,disp_X,disp_Z,disp_mesh,writefig,runFolder)
%% this function plot values in/of domain

%%  input
%%  domain
%%  disp_edensity,disp_coll_freq,disp_mu,disp_kappa,disp_X,disp_Z,disp_mesh
%%  output
%%  none
%%  internal parameters
nb_lines=10;

cd (runFolder);

%% larger LOS line
%     for it=1:3
%         scdir_x(it)=pooo1(1)+cos(scdir/180*pi)*(it-1)*5;                 scdir_y(it)=pooo1(2)+sin(scdir/180*pi)*(it-1)*5;
%         scdirr1_x(it)=pooo1(1)+cos((scdir+scdir_range)/180*pi)*(it-1)*3; scdirr1_y(it)=pooo1(2)+sin((scdir+scdir_range)/180*pi)*(it-1)*3;
%         scdirr2_x(it)=pooo1(1)+cos((scdir-scdir_range)/180*pi)*(it-1)*3;	scdirr2_y(it)=pooo1(2)+sin((scdir-scdir_range)/180*pi)*(it-1)*3;
%     end


if disp_edensity  %% plot electrondensity/critical density
    figure
    hold on
    
    VD=[criticaldensity criticaldensity];
    
    
    for z=1:domain.nozones
        % plot boundary
        plot( domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:)) , ...
            domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:)) , ...
            'k-','LineWidth',2)
        
        
        % plot electrondensity countours
        %                [C,h]=tricontour(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  ),nb_lines);
        %[C,h]=tricontour(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  ),VE);
        %[C,h]=tricontour(    domain.( strcat('zone',num2str(1)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(1)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(1)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(1)) ).variables(domain.nova-5,:)  ),VE);
        %VE=[1e13,1e14,1e15,1e16,1e17,1e18];
        
        
        h=trisurf( domain.( strcat('zone',num2str(z))).delaunay     , ...
            squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  , ...
            squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  , ...
            squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  ) );
        
        xlabel('X [m]');ylabel('Y [m]');zlabel('Z');
        set(h, 'edgecolor','none');
        set(gcf, 'Position', [10, 10, 1100, 800]);
        
        
        %   xlim([-0.17 6]);ylim([0 5.185]);caxis([1e16 1e17]);
        maxX(z) = max(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        minX(z) = min(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        maxY(z) = max(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        minY(z) = min(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        maxC(z) = max(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  );
        minC(z) = min(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  );
        maxX = max(maxX);
        minX = min(minX);
        maxY = max(maxY);
        minY = min(minY);
        maxC = max(maxC);
        minC = min(minC);
        axis equal; xlim([minX maxX]); ylim ([minY maxY]); caxis([minC maxC]);
        
        
        %h.Fill='on'
        % %plot critical electrondensity, Davis page 61, eq. 2.54 and k of page 71
        %{
                if max_edensity > criticaldensity
                    [D,i]=tricontour(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-5,:)  ));
                    i(1).LineWidth=2;
                    i(1).EdgeColor='k';
                    %text(0,5,strcat('critical density=',num2str(criticaldensity,'%1.1e')));
                end
        %}
        
        hcb=colorbar;title(hcb,'X_e [1/m^3]          ');
        box on;grid off;
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
        
    end
    hold off
    
    if writefig
        savefig('e_density')
        %print('-painters','e_density','-depsc')
        %print('-painters','e_density','-dpng','-r300')
        print('-painters','e_density','-dpng')
        
        %ax = gca;
        % Requires R2020a or later
        %exportgraphics(ax,'e_density.png','Resolution',200)
    end
    
end

if disp_coll_freq    %%plot collision surface plot
    figure
    hold on
    
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
        
        h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)  ) );
        set(h, 'edgecolor','none');
        set(gcf, 'Position', [10, 10, 1100, 800]);
        maxX(z) = max(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        minX(z) = min(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        maxY(z) = max(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        minY(z) = min(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        maxC(z) = max(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)  );
        minC(z) = min(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-4,:)  );
        maxX = max(maxX);
        minX = min(minX);
        maxY = max(maxY);
        minY = min(minY);
        maxC = max(maxC);
        minC = min(minC);
%        axis equal; xlim([minX maxX]); ylim ([minY maxY]); caxis([minC maxC]);
        
        hcb=colorbar;title(hcb,'collision freq. [1/s]                 ');
        box on;grid off;
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
        
    end
    hold off
    
    if writefig
        savefig('coll_freq')
        %print('-painters','coll_freq','-depsc')
        %print('-painters','coll_freq','-dpng','-r300')
        print('-painters','coll_freq','-dpng')
        %ax = gca;
        % Requires R2020a or later
        %exportgraphics(ax,'coll_freq.png','Resolution',200)
    end
    
end

if disp_X    %% plot X
    figure
    hold on
    VX0=[1 10^7];
    for z=1:domain.nozones
        % plot boundary
        plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)
        % plot X
        %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)),nb_lines);
        %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)),VX);
        %VX=[1e-17,1e-15,1e-13,1e-11,1e-10,1e-8,1e-6,1e-4,1e-2,1e-1,2e-1:0.02:1-2e-1,1-1e-1,1-1e-3,1,3,5,7,9,20,40,60,80,100];
        
        h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)  ) );
        set(h, 'edgecolor','none');
        set(gcf, 'Position', [10, 10, 1100, 800]);
        maxX(z) = max(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        minX(z) = min(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        maxY(z) = max(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        minY(z) = min(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        maxC(z) = max(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)  );
        minC(z) = min(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)  );
        maxX = max(maxX);
        minX = min(minX);
        maxY = max(maxY);
        minY = min(minY);
        maxC = max(maxC);
        minC = min(minC);
%        axis equal; xlim([minX maxX]); ylim ([minY maxY]); caxis([minC maxC]);
        %xlim([-0.35 6.1]);ylim([0 5.2]);caxis([0 10])
        
        hcb=colorbar;title(hcb,'X [-]');
        box on;grid off;
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
        if max_edensity > criticaldensity
            [D,i]=tricontour(    domain.( strcat('zone',num2str(z)) ).delaunay    ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-3,:)  ),VX0);
            i(1).LineWidth=2;
            i(1).EdgeColor='k';
            %text(0,4.5,'X=1');
        end
        
    end
    hold off
    
    if writefig
        savefig('X')
        %print('-painters','X','-depsc')
        %print('-painters','X','-dpng','-r300')
        print('-painters','coll_freq','-dpng')
        %ax = gca;
        % Requires R2020a or later
        %exportgraphics(ax,'X.png','Resolution',200)
    end
    
end

if disp_Z    %% plot Z
    figure
    %axis equal;
    hold on
    
    for z=1:domain.nozones
        %% plot boundary
        plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)
        %% plot Z
        %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-2,:)),nb_lines);
        %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay  ,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-2,:)),VZ);
        %VZ=[1e-17,1e-16,1e-15,1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6];
        
        h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-2,:)  ) );
        set(h, 'edgecolor','none');
        set(gcf, 'Position', [10, 10, 1100, 800]);
        maxX(z) = max(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        minX(z) = min(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        maxY(z) = max(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        minY(z) = min(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        maxC(z) = max(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-2,:)  );
        minC(z) = min(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-2,:)  );
        maxX = max(maxX);
        minX = min(minX);
        maxY = max(maxY);
        minY = min(minY);
        maxC = max(maxC);
        minC = min(minC);
%        axis equal; xlim([minX maxX]); ylim ([minY maxY]); caxis([minC maxC]);
        %xlim([-0.35 6.1]);ylim([0 5.2]);caxis([0 1e-6])
        
        
        hcb=colorbar;title(hcb,'Z [-]       ');
        box on;grid off;
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
    end
    hold off
    
    if writefig
        
        savefig('Z')
        %print('-painters','Z','-dpng','-r300')
        print('-painters','Z','-dpng')
        %print('-painters','Z','-depsc')
        
        %ax = gca;
        % Requires R2020a or later
        %exportgraphics(ax,'Z.png','Resolution',200)
    end
    
end

if disp_mu    %% plot refractive index
    figure
    hold on
    for z=1:domain.nozones
        % % plot boundary
        plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)
        % % plot ri
        %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),nb_lines);
        %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),VRI);
        VRI=[1e-17,1e-15,1e-13,1e-11,1e-9,1e-7,1e-5,1e-3,1e-1,2e-1:0.02:1-2e-2,1-1e-1,1-1e-3,1-1e-5,1-1e-7,1-1e-9,1-1e-11,1-1e-13,1-1e-15,1-1e-17,1];
        
        h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)  ) );
        set(h, 'edgecolor','none');
        set(gcf, 'Position', [10, 10, 1100, 800]);
        %xlim([-0.35 6.1]);ylim([0 5.2]);caxis([0 1])
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
%        axis equal; xlim([minX maxX]); ylim ([minY maxY]); caxis([minC maxC]);
        
        hcb=colorbar;title(hcb,'\mu [-]');
        box on;grid off;
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
    end
    %img=imread('cases/ExoMars/ExoMarsHalfCapsule.png');
    %imagesc([1 1],[5 5],img);
    hold off
    
    if writefig
        savefig('mu')
        %print('-painters','mu','-dpng','-r300')
        print('-painters','mu','-dpng')
        
        %print('-painters','mu','-depsc')
        
        %ax = gca;
        % Requires R2020a or later
        %exportgraphics(ax,'mu.png','Resolution',200)
        %print('-painters','mu','-depsc')
        %imwrite(fh,'mu_s.png')
    end
    
end

if disp_kappa    %% plot absorption coefficient
    
    fh = figure();
    % fh.WindowState = 'maximized'; %maximize window
    hold on
    for z=1:domain.nozones
        % plot boundary
        plot(      domain.( strcat('zone',num2str(z)) ).variables(1,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  domain.( strcat('zone',num2str(z)) ).variables(2,domain.( strcat('zone',num2str(z)) ).bound(:))  ,  'k-','LineWidth',2)
        % plot kappa
        %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),nb_lines);
        %[C,h]=tricontour(domain.( strcat('zone',num2str(z)) ).delaunay,squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:)),squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova-1,:)),VRI);
        VRI=[1e-17,1e-15,1e-13,1e-11,1e-9,1e-7,1e-5,1e-3,1e-1,2e-1:0.02:1-2e-2,1-1e-1,1-1e-3,1-1e-5,1-1e-7,1-1e-9,1-1e-11,1-1e-13,1-1e-15,1-1e-17,1];
        
        h=trisurf(            domain.( strcat('zone',num2str(z))).delaunay     ,   squeeze(domain.( strcat('zone',num2str(z)) ).variables(1,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(2,:))  ,  squeeze(domain.( strcat('zone',num2str(z)) ).variables(domain.nova,:)  ) );
        set(h, 'edgecolor','none');
        set(gcf, 'Position', [10, 10, 1100, 800]);
        maxX(z) = max(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        minX(z) = min(domain.( strcat('zone',num2str(z)) ).variables(1,:));
        maxY(z) = max(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        minY(z) = min(domain.( strcat('zone',num2str(z)) ).variables(2,:));
        maxC(z) = max(domain.( strcat('zone',num2str(z)) ).variables(domain.nova,:)  );
        minC(z) = min(domain.( strcat('zone',num2str(z)) ).variables(domain.nova,:)  );
        maxX = max(maxX);
        minX = min(minX);
        maxY = max(maxY);
        minY = min(minY);
        maxC = max(maxC);
        minC = min(minC);
       % axis equal; xlim([minX maxX]); ylim ([minY maxY]); caxis([minC maxC]);
        
        hcb=colorbar;title(hcb,'\kappa [1/m]');
        box on;grid off;
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
    end
    %img=imread('cases/ExoMars/ExoMarsHalfCapsule.png');
    %imagesc([1 1],[5 5],img);
    hold off
    
    if writefig
        savefig('kappa')
        %print('-painters','kappa','-dpng','-r300')
        print('-painters','kappa','-dpng','-r300')
        
        %print('-painters','kappa','-depsc')
        
        %ax = gca;
        % Requires R2020a or later
        %exportgraphics(ax,'kappa.png','Resolution',200)
    end
    
end

if disp_mesh    %%plot mesh
    
    figure
    axis equal;
    hold on
    for z=1:domain.nozones
        triplot( domain.( strcat('zone',num2str(z)) ).delaunay  ,  domain.( strcat('zone',num2str(z)) ).variables(1,:)'   ,   domain.( strcat('zone',num2str(z)) ).variables(2,:)'   );
    end
    hold off
    if writefig
        savefig('mesh')
        %print('-painters','mesh','-depsc')
        %print('-painters','mesh','-dpng','-r300')
        print('-painters','mesh','-dpng')
        
    end
end

cd ..;

end % function plotdomain


function [c,h]=tricontour(tri,x,y,z,nv)
%TRICONTOUR Triangular Contour Plot.
% TRICONTOUR(TRI,X,Y,Z,N) draws scalar N contour lines treating the values
% in Z as heights above a plane. TRI,X,Y,and Z define a triangulation where
% the triangles are defined by the M-by-3 face matrix TRI, such as that
% returned by DELAUNAY. Each row of TRI contains indices into the X,Y, and
% Z vertex vectors to define a single triangular face. Contours are
% computed directly from the triangulation rather than interpolating back
% to a cartesian grid using GRIDDATA.
% TRICONTOUR(TRI,X,Y,Z,V) draws length(V) contour lines at the values
% specified in vector V.
% TRICONTOUR(TRI,X,Y,Z,[v v]) draws a single contour line at the level v.
%
% [C,H] = TRICONTOUR(...) returns contour matrix C as described in CONTOURC
% and a vector of handles H to the created patch objects.
% H can be used to set patch properties.
% CLABEL(C) or CLABEL(C,H) labels the contour levels.
%
% Example:
%           x=linspace(-3,3,39);
%           y=linspace(-2.5,2.5,49);
%           [xx,yy]=meshgrid(x,y);
%           zz=peaks(xx,yy);
%           v=-3:1:5; % contour levels
%           subplot(1,2,1)
%           [C,h]=contour(xx,yy,zz,v);   % standard contour for comparison
%           clabel(C)
%           title Contour
%
%           idx=randperm(numel(zz));     % grab some scattered indices
%           n=idx(1:ceil(numel(zz)/2))'; % one half of them
%           x=xx(n);                     % get scattered data
%           y=yy(n);
%           z=zz(n);
%           tri=delaunay(x,y);           % triangulate scattered data
%           subplot(1,2,2)
%           [C,h]=tricontour(tri,x,y,z,v);
%           clabel(C,h)
%           title TriContour
%
% view(3) displays the contour in 3-D.
%
% See also DELAUNAY, CONTOUR, TRIMESH, TRISURF, TRIPLOT, PATCH.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2006-05-07, 2006-05-16, 2006-07-25
d = nargin;
if nargin<5
    error('Not Enough Input Arguments.')
end
x=x(:);	% convert input data into column vectors
y=y(:);
z=z(:);
xlen=length(x);
if ~isequal(xlen,length(y),length(z))
    error('X, Y, and Z Must Have the Same Number of Elements.')
end
if size(tri,2)~=3 || any(tri(:)<0) || any(tri(:)>xlen)
    error('TRI Must Be a Valid Triangulation of the Data in X, Y, Z.')
end

zs=z(tri);
zmax=max(max(zs));              % find max and min in z data that is in tri
zmin=min(min(zs));

if length(nv)==1                                 % nv is number of contours
    zlev=linspace(zmax,zmin,nv+2);
elseif length(nv)==2 && nv(1)==nv(2)              % nv is one contour level
    zlev=nv(1);
else                                       % nv is vector of contour levels
    zlev=sort(nv,'descend');
end
zlev(zlev>=zmax | zlev<=zmin)=[];  % eliminate contours outside data limits
nlev=length(zlev);

if nlev==0
    error('No Contours to Plot. Chosen Contours Outside Limits of Data.')
end

% precondition the input data
[zs,zidx]=sort(zs,2);         % sort vertices by z value ascending
for k=1:size(zs,1)            % shuffle triangles to match
    tri(k,:)=tri(k,zidx(k,:));
end

hax=newplot;                  % create new axis if needed
h=[];                         % patch handle storage
C=zeros(2,0);                 % Clabel data storage
cs=[2 1];                     % column swap vector cs(1)=2, cs(2)=1;

% Main Loop ---------------------------------------------------------------
for v=1:nlev                  % one contour level at a time
    zc=zlev(v);                % chosen level
    above=zs>=zc;              % true for vertices above given contour
    numabove=sum(above,2);     % number of triangle vertices above contour
    tri1=tri(numabove==1,:);   % triangles with one vertex above contour
    tri2=tri(numabove==2,:);   % triangles with two vertices above contour
    n1=size(tri1,1);           % number with one vertex above
    n2=size(tri2,1);           % number with two vertices above
    
    edge=[tri1(:,[1 3])        % first column is indices below contour level
        tri1(:,[2 3])        % second column is indices above contour level
        tri2(:,[1 2])
        tri2(:,[1 3])];
    if n1==0                   % assign edges to triangle number
        n=[1:n2 1:n2]';
    elseif n2==0
        n=[1:n1 1:n1]';
    else
        n=[1:n1 1:n1 n1+(1:n2) n1+(1:n2)]';
    end
    
    [edge,idx]=sortrows(edge);    % put shared edges next to each other
    n=n(idx);                     % shuffle triangle numbers to match
    
    idx=all(diff(edge)==0,2);     % find shared edges
    idx=[idx;false]|[false;idx];  % True for all shared edges
    
    % eliminate redundant edges, two triangles per interior edge
    edgeh=edge(~idx,:);           % hull edges
    nh=n(~idx);                   % hull triangle numbers
    if ~isempty(nh)
        nh(end,2)=0;               % zero second column for hull edges
    end
    edges=edge(idx,:);            % shared edges
    edges=edges(1:2:end-1,:);     % take only unique edges
    ns=n(idx);                    % interior triangle numbers
    ns=[ns(1:2:end) ns(2:2:end)]; % second column is second triangle
    edge=[edgeh;edges];           % unique edges
    nn=[nh;ns];                   % two columns of triangle numbers
    ne=size(edge,1);              % number of edges
    
    flag=true(ne,2);              % true for each unused edge per triangle
    tmp=zeros(ne+1,1);            % contour data temporary storage
    
    xe=x(edge);                   % x values at vertices of edges
    ye=y(edge);                   % y values at  vertices of edges
    ze=z(edge);                   % z data at  vertices of edges
    
    alpha=(zc-ze(:,1))./(ze(:,2)-ze(:,1)); % interpolate all edges
    xc=alpha.*(xe(:,2)-xe(:,1)) + xe(:,1); % x values on this contour
    yc=alpha.*(ye(:,2)-ye(:,1)) + ye(:,1); % y values on this contour
    
    while any(flag)	% while there are still unused edges -----------------
        
        xtmp=tmp;
        ytmp=tmp;
        [ir,ic]=find(flag,1);            % find next unused edge
        flag(ir,ic)=false;               % mark this edge used
        
        k=1;                             % first data point in subcontour
        xtmp(k)=xc(ir);                  % store data from this edge
        ytmp(k)=yc(ir);
        
        while true     % complete this subcontour ---------------------------
            
            [ir,ic]=find(flag&nn(ir,ic)==nn,1);% find other edge of triangle
            flag(ir,ic)=false;            % mark this edge used
            k=k+1;
            xtmp(k)=xc(ir);               % store data from this edge
            ytmp(k)=yc(ir);
            
            ic=cs(ic);                    % other triangle that shares edge
            
            if nn(ir,ic)==0               % reached hull, subcontour complete
                k=k+1;
                xtmp(k)=nan;               % don't let subcontour close
                ytmp(k)=nan;
                break
            elseif ~flag(ir,ic)           % complete closed subcontour
                break
            else                          % more points remain on subcontour
                flag(ir,ic)=false;         % mark this edge used
            end
        end % while true ----------------------------------------------------
        xtmp(k+1:end)=[];                % throw away unused storage
        ytmp(k+1:end)=[];                % xtmp,ytmp contain subcontour
        
        if nargout<2                     % plot the subcontour
            patch('XData',xtmp,'YData',ytmp,'CData',repmat(zc,k,1),...
                'Parent',hax,'FaceColor','none','EdgeColor','flat',...
                'UserData',zc)
            C=horzcat(C,[zc xtmp';k ytmp']); % contour label data
        else                             % plot subcontour and create output
            h=[h;patch('XData',xtmp,'YData',ytmp,'CData',repmat(zc,k,1),...
                'Parent',hax,'FaceColor','none','EdgeColor','flat',...
                'UserData',zc)]; %#ok
            C=horzcat(C,[zc xtmp';k ytmp']); % contour label data
        end
    end % while any(flag) --------------------------------------------------
end % for v=1:nlev
if nargout
    c=C;
end
end
