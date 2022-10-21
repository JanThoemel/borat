tecplot_home = '/Applications/Tecplot 360 EX 2020 R2/Tecplot 360 EX 2020 R2.app/Contents';
tecio_path = strcat(tecplot_home, '/Frameworks/libtecio.dylib');
tecio_header_path = '/TecIO/TECIO.h'; %Uses simplified TECIO.h included in package.



if ~libisloaded('tecio')
    [notfound, warnings] = loadlibrary(tecio_path, tecio_header_path,...
        'alias', 'tecio');
end


%%%%%%%%%%% Output data to .szplt data file. %%%%%%%%%%%
output_fname = 'Rays.szplt';

dataset_title = 'raytracing 3D ';
var_names = 'x0, x1, x2, phase, path, n';
file_format = 1;  % 0 - .plt;  1 - .szplt
file_type   = 0;  % 0 - grid & solution; 1 - grid only; 2 - solution only
data_type   = 2;  % 1 - single; 2 - double; ...

total_zones=size(rays_solution_3D(:,1,1),1);

zone_title = 'Ray n ';

[isok,~,~,~,~,filehandle] = calllib('tecio', 'tecFileWriterOpen', ...
    output_fname, dataset_title, var_names, ...
    file_format, file_type, data_type, [],[] );

for nzone=1:total_zones
    
    total_points=nnz(rays_solution_3D(nzone,:,1));
    I=total_points;
    J=1;
    K=1;
    zname = [ zone_title num2str(nzone,'%02d')];
    
    [isok,~,~,~,~,~,~,n_zone] = calllib('tecio', 'tecZoneCreateIJK', ...
        filehandle, zname, I, J, K, [], [], [],[],0,0,0,nzone);
    
     vars= [rays_solution_3D(nzone,1:total_points,1)' rays_solution_3D(nzone,1:total_points,2)' rays_solution_3D(nzone,1:total_points,3)' ...
         rays_solution_3D(nzone,1:total_points,7)' rays_solution_3D(nzone,1:total_points,8)' rays_refractive_index(nzone,1:total_points,1)'];

    for v = 1:size(vars,2)
        isok = calllib('tecio', 'tecZoneVarWriteDoubleValues', ...
            filehandle, nzone, v, 1, total_points, vars(:,v));
    end
    
end %nzone

calllib('tecio', 'tecFileWriterClose', filehandle);

if libisloaded('tecio')
    unloadlibrary('tecio')
end

