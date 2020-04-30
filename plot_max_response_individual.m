function plot_max_response_individual(binnedC, comparison, neuron2plot,t)
%% plot bar chart for maximum response bin
temp1 = squeeze(nanmean(binnedC(:,comparison,:,neuron2plot))); %find mean binnedC for each bin
temp2 = squeeze(nanstd(binnedC(:,comparison,:,neuron2plot)));
temp3 = squeeze(binnedC(:,comparison,:,neuron2plot));
for i =1:size(binnedC, 3)
    AA(i) = length(find(~isnan(temp3(:,i))));
end
sem = temp2./sqrt(AA);
% t=categorical({'100', '75', '55', '45', '25', '0'});



    for i=1:6
        hold on
        bar(t(i), temp1(i));
        errorbar(t(i),temp1(i),sem(i));
        
    end
    title(append('Neuron ',num2str(neuron2plot)));
    ylabel('Inferred Spiking (AU)')
    hold off
   
end