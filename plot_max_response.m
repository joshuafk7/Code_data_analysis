function plot_max_response(binnedC, comparison, sig_neurons,t)
%% plot bar chart for maximum response bin
figure
temp1 = squeeze(nanmean(binnedC(:,comparison,:,sig_neurons))); %find mean binnedC for each bin
% temp2 = squeeze(nanstd(binnedC(:,3,:,taste_resp_neurons)));

% t=categorical({'100', '75', '55', '45', '25', '0'});
sz = round(sqrt(length(sig_neurons)))+1; %pick number of subplots

for j=1:length(sig_neurons)
    subplot(sz,sz-1,j)
    for i=1:length(t)
        hold on
        bar1 = bar(t(i), temp1(i,j));
%         errorbar(t(i),temp1(i,j),temp2(i,j));
        title(append('Neuron ',num2str(sig_neurons(j))));
        ylabel('Inferred Spiking')
%         xlabel('Taste')
    end
    hold off
   
end
end