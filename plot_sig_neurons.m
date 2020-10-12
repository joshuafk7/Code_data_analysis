function plot_sig_neurons(sig_neurons, plotbins, plotbinnedC2,summary)
%plot trial averages for responsive neurons 
%inputs: sig_neurons - a list of significantly modulated neurons
%plotbins: bin sizes for plotting
%plotbinnedC2: binned data for plotting
%summary:contains behavior info including animal ID, date, etc...
figure(7539)
sz = round(sqrt(length(sig_neurons)))+1; %pick number of subplots

for j=1:length(sig_neurons)
    subplot(sz,sz-1,j) %plot all of the taste responsive neurons
    plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,sig_neurons(j)))) %average across trials and tastes
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(sig_neurons(j))))
end
suptitle(summary.mouseID);
end