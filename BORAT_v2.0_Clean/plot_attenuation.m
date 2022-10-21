function plot_attenuation(final_dir_attenuation,writefig,runFolder)
    
    cd (runFolder);

    final_dir_attenuation=sortrows(final_dir_attenuation);
    
    
    
    figure
    polarplot( [final_dir_attenuation(:,1)'*pi/180 final_dir_attenuation(1,1)*pi/180], [final_dir_attenuation(:,2)' final_dir_attenuation(1,2)] );
    hold on
    hold off
    if writefig
        savefig('polar_mu')
        print('-painters','polar_mu','-depsc')
        print('-painters','polar_mu','-dpng')
    end
    
    
    cd ..;
    
end