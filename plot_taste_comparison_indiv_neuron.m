function plot_taste_comparison_indiv_neuron(neuron2plot,plotbins, plotbinnedC)
    subplot(3,1,1)
    hold on
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,1,neuron2plot,1)))
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,6,neuron2plot,1)))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(neuron2plot)))  
    ylabel('Inferred Spikes (AU)')
    legend('Suc','NaCl','Location','Northwest')
    hold off
    
    subplot(3,1,2)
    hold on
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,2,neuron2plot,1)))
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,5,neuron2plot,1)))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(neuron2plot)))
    ylabel('Inferred Spikes (AU)')
    legend('75','25','Location','Northwest') 
    hold off

    
  subplot(3,1,3)
  hold on
  plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,3,neuron2plot,1)))
  plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,4,neuron2plot,1)))
  y=ylim;
  line([0 0],[0 y(2)],'Color','r','LineStyle','--')
  title(append('Neuron ',num2str(neuron2plot)))
  legend('55','45','Location','Northwest') 
    
    
    title(append('Neuron ',num2str(neuron2plot)))
    xlabel('Time (sec)')
    ylabel('Inferred Spikes (AU)')
    hold off
end