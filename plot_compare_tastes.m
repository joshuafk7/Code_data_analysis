function plot_compare_tastes(plotbins, plotbinnedC, CC)
titles = {'Suc vs NaCl','S\_75\_25 vs S\_25\_75','S\_55\_45 vs S\_45\_55'};
for i=1:3
    figure
    suptitle(titles(i))
    ZZ = find(CC(:,2)==i);
   sz = round(sqrt(length(ZZ)))+1; %pick number of subplots
for j=1:length(ZZ)
    subplot(sz,sz-1,j)
    hold on
    if i==1
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,1,CC(ZZ(j),1)),1))
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,6,CC(ZZ(j),1)),1))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(CC(ZZ(j),1))))

    if j==1 
        legend('Suc','NaCl','Location','Northwest')
    end
    end
    if i==2
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,2,CC(ZZ(j),1)),1))
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,5,CC(ZZ(j),1)),1))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(CC(ZZ(j),1))))
    if j==1 
        legend('75','25','Location','Northwest') 
    end

    end
    if i==3
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,3,CC(ZZ(j),1)),1))
    plot(plotbins(1:end-1), nanmean(plotbinnedC(:,:,4,CC(ZZ(j),1)),1))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(CC(j,1))))
    if j==1
        legend('55','45','Location','Northwest') 
    end
    end
    title(append('Neuron ',num2str(CC(ZZ(j),1))))
    hold off
end

end
end