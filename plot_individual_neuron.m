function ylimit = plot_individual_neuron(n,neuron2plot, plotbins, plotbinnedC2)
    switch n
        case 1
    plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,neuron2plot))) %average across trials and tastes
    y=ylim;
    ylimit = y(2);
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(neuron2plot)))
    xlabel('Time (sec)')
    ylabel('Inferred Spiking (AU)')
        case 2
            hold on
           plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,1,neuron2plot)))
            plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,2,neuron2plot)))%average across trials and tastes
    y=ylim;
    ylimit = y(2);
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    rectangle('Position',[4 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5]);
    title(append('Neuron ',num2str(neuron2plot)));
    xlabel('Time (sec)')
    ylabel('Inferred Spiking (AU)') 
            
    end
end