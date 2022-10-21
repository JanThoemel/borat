function surfaceangle=isolinefinding(domain,snell_po_base,ss)

%% 	find electron iso surface 
%%
%% input
%%
%% output
%%
%% remarks
%! experiment with different/adaptive circleradius and number of perimeter angles
%! do something about 1) totally flat surfaces 2) strange surfaces (,e.g. interpolation)
%% internal parameters
    circleradius=ss/7;
    numberofcircledivisions=360;

    %% draw circle with radius ss around po_snellbase and determine ri
    perimeter=zeros(numberofcircledivisions,4);
    perimeter(:,1)=1:360/numberofcircledivisions:numberofcircledivisions;    %angles
    
    ri_snellbase=interpolation(domain,snell_po_base,domain.nova-1);

    for a=1:size(perimeter,1)
	   perimeter(a,2)=snell_po_base(1)+cos(perimeter(a,1)/180*pi)*circleradius; %x
	   perimeter(a,3)=snell_po_base(2)+sin(perimeter(a,1)/180*pi)*circleradius; %y
	   perimeter(a,4)=interpolation(domain,perimeter(a,2:3),domain.nova-1); %ri

	   %this is sometimes discontinous, likely poor interpolation on a poor mesh effect
	   %!check it        
    end
    
    %!check uniformity amd continuity for ri(over circle)            
    %for debugging:
	  % plot(perimeter(:,1),perimeter(:,4));
	  % pause
        
    %%find on perimeter the two points that are closest to ri i.e. find
    %%indizes of perimeter, this defines the ri-isosurface
    a1=1;a2=360;    
    diffri1 =abs(perimeter(1,4)-ri_snellbase);
    for a=1:1:179
		if abs(perimeter(a,4)-ri_snellbase)<diffri1
			diffri1=abs(perimeter(a,4)-ri_snellbase);
			a1=a;
		end
    end
 
    diffri2 =abs(perimeter(360,4)-ri_snellbase);
    for a=359:-1:180
		if abs(perimeter(a,4)-ri_snellbase)<diffri2
			diffri2=abs(perimeter(a,4)-ri_snellbase);
			a2=a;
		 end
    end
    surfaceangle=a1;

%{
    %as of 25th sept 2017 all warnings/errors have been verified not to impact, therefore the
    %following check is disabled.
    %% check flatness of surface, warn or terminate depending on level non-flatness
    %% if abs(a2-a1)>356 then the determination of the isosurface fails
    %% likely because ri is constant over this circle 
%    if (abs(a2-a1)<170 || abs(a2-a1)>190) && abs(a2-a1)<356 
        p = polyfit(perimeter(:,1),perimeter(:,4),4)
        %[p,S,mu] = polyfit(x,y,n) also returns mu, which is a two-element vector with centering and scaling values. mu(1) is mean(x), and mu(2) is std(x). Using these values, polyfit centers x at zero and scales it to have unit standard deviation
        fitter=polyval(p,perimeter(:,1));
        wendep1=0.25*p(2)/p(1)+sqrt((0.25*p(2)/p(1))^2-1/6*p(3)/p(1));
        wendep2=0.25*p(2)/p(1)-sqrt((0.25*p(2)/p(1))^2-1/6*p(3)/p(1));
        plot(perimeter(:,1),perimeter(:,4)); %plot ri over angles
        hold on
        t=text(0,ri_snellbase,strcat('a1=',num2str(a1),' a2=',num2str(a2))); %show in figure value of two angles
        plot(perimeter(:,1),fitter);
        warning('snellslaw - ri isosurface not well defined, a1: %.2f, a2: %.2f, a2-a1:%.2f',a1,a2,a2-a1)
        pause
        close;
 %   end %if        
%}
    %{
    if (abs(a2-a1)<140 || abs(a2-a1)>220) && abs(a2-a1)<356 
        %plot(perimeter(:,1),perimeter(:,4)); %plot ri over angles
        %t=text(0,ri_snellbase,strcat('a1=',num2str(a1),' a2=',num2str(a2))); %show in figure value of two angles
		%error('snellslaw - ri isosurface insuffciently defined, a1: %.2f, a2: %.2f, a2-a1:%.2f',a1,a2,a2-a1)
    end %if        
%}
    
end %isoline finding
