clc; clear all; close all;  initime=cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% emission of 1 angle, comparison between eikonal and snell law
load('CFDsolutions/raytracing_solution_62.mat')

runFolder='Solutions';

dir1=30; %angle of the ray
%% Stepsize

ss=0.01;
ss1=0.1;
ss2=0.05;
ss3=0.02;
ss4=0.01;

%% initialization

dir2=dir1;
maxangles=1;

itpo1=itpo;
itpo2=itpo;
itpo3=itpo;
itpo4=itpo;

%% plot domain

%% plot solution
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
    %xlim([-0.35 6.1]);ylim([0 5.2]);caxis([0 1])
    axis equal; xlim([-0.5 16]); ylim ([-6 6]); caxis([0 1])
    axis tight;


    hcb=colorbar;title(hcb,'\mu [-]');
    box on;grid off;
    set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
end
%img=imread('cases/ExoMars/ExoMarsHalfCapsule.png');
%imagesc([1 1],[5 5],img);
hold on
colormap('gray');

%% eikonal integration

tic
[itdir, itpo,symmetrylineencounter]=eikonal2D(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, absorptionlimits(2)-10,symmetryline,ss);
fprintf('\n  ');
toc

z_max = max(max(get(h,'Zdata')));
scdir_z=[z_max,z_max,z_max];

%%plot rays in domain
itpoz=z_max*ones(   size( itpo(1,:,1),2)   ,1);%z coordinates of lines to be plotted
for a=1:maxangles
    %!can they be colored according to attenuation?
    plot3(  itpo(1,1:nnz(itpo(1,:,a)),a)  ,  itpo(2,1:nnz(itpo(1,:,a)),a),itpoz(1:nnz(itpo(1,:,a))),'bo','LineWidth',1);
end

hold on

%% case1 Snell
iss=0.1;

tic
[itdir2, itpo4,symmetrylineencounter]=raytracing_snell_law(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, 100,symmetryline,ss1);
fprintf('\n  ');
toc

z_max = max(max(get(h,'Zdata')));
scdir_z=[z_max,z_max,z_max];

%%plot rays in domain
itpoz=z_max*ones(   size( itpo4(1,:,1),2)   ,1);%z coordinates of lines to be plotted
for a=1:maxangles
    %!can they be colored according to attenuation?
    plot3(  itpo4(1,1:nnz(itpo4(1,:,a)),a)  ,  itpo4(2,1:nnz(itpo4(1,:,a)),a),itpoz(1:nnz(itpo4(1,:,a))),'ro','LineWidth',1);
end

hold on


%% case2 Snell
tic
[itdir2, itpo3,symmetrylineencounter]=raytracing_snell_law(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, 100,symmetryline,ss2);
fprintf('\n  ');
toc

z_max = max(max(get(h,'Zdata')));
scdir_z=[z_max,z_max,z_max];

%%plot rays in domain
itpoz=z_max*ones(   size( itpo3(1,:,1),2)   ,1);%z coordinates of lines to be plotted
for a=1:maxangles
    %!can they be colored according to attenuation?
    plot3(  itpo3(1,1:nnz(itpo3(1,:,a)),a)  ,  itpo3(2,1:nnz(itpo3(1,:,a)),a),itpoz(1:nnz(itpo3(1,:,a))),'go','LineWidth',1);
end

hold on

%% case3 Snell

tic
[itdir1, itpo1,symmetrylineencounter]=raytracing_snell_law(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, 100,symmetryline,ss3);
fprintf('\n  ');
toc

z_max = max(max(get(h,'Zdata')));
scdir_z=[z_max,z_max,z_max];

%%plot rays in domain
itpoz=z_max*ones(   size( itpo1(1,:,1),2)   ,1);%z coordinates of lines to be plotted
for a=1:maxangles
    %!can they be colored according to attenuation?
    plot3(  itpo1(1,1:nnz(itpo1(1,:,a)),a)  ,  itpo1(2,1:nnz(itpo1(1,:,a)),a),itpoz(1:nnz(itpo1(1,:,a))),'ko','LineWidth',1);
end

hold on

%% case4 Snell

tic
[itdir2, itpo2,symmetrylineencounter]=raytracing_snell_law(domain,pooo1,pooo2,dir1,dir2,maxsteps,maxangles, 100,symmetryline,ss4);
fprintf('\n  ');
toc

z_max = max(max(get(h,'Zdata')));
scdir_z=[z_max,z_max,z_max];

%%plot rays in domain
itpoz=z_max*ones(   size( itpo2(1,:,1),2)   ,1);%z coordinates of lines to be plotted
for a=1:maxangles
    %!can they be colored according to attenuation?
    plot3(  itpo2(1,1:nnz(itpo2(1,:,a)),a)  ,  itpo2(2,1:nnz(itpo2(1,:,a)),a),itpoz(1:nnz(itpo2(1,:,a))),'mo','LineWidth',1);
end

hold off


%%
function plotdomainmu(domain)

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
        %xlim([-0.35 6.1]);ylim([0 5.2]);caxis([0 1])
        axis equal; xlim([-0.5 16]); ylim ([-6 6]); caxis([0 1])
        axis tight;


        hcb=colorbar;title(hcb,'\mu [-]');
        box on;grid off;
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontName','Times New Roman', 'Fontweight','bold','Linewidth',2,'FontSize',12,'TickLength',[0.02, 0.002])
    end
    %img=imread('cases/ExoMars/ExoMarsHalfCapsule.png');
    %imagesc([1 1],[5 5],img);
    hold off
    colormap('gray');

end % function plotdomain

function [c,h]=tricontour(tri,x,y,z,nv)


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
