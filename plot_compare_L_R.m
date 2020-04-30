function plot_compare_L_R(plotbins, plotbinnedC3, CC)
titles = {'Correct','Error'};
for i=1:2
    figure
    suptitle(titles(i));
    ZZ = find(CC(:,2)==i);
   sz = round(sqrt(length(ZZ)))+1; %pick number of subplots
for j=1:length(ZZ)
    subplot(sz,sz-1,j)
    hold on
    if i==1
        
    plot(plotbins(1:end-1), nanmean(plotbinnedC3(:,:,1,CC(ZZ(j),1)),1))
    plot(plotbins(1:end-1), nanmean(plotbinnedC3(:,:,2,CC(ZZ(j),1)),1))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
    title(append('Neuron ',num2str(CC(ZZ(j),1))))

    if j==1 
        legend('Left','Right','Location','Northwest')
    end
    end
    if i==2
        
    plot(plotbins(1:end-1), nanmean(plotbinnedC3(:,:,3,CC(ZZ(j),1)),1))
    plot(plotbins(1:end-1), nanmean(plotbinnedC3(:,:,4,CC(ZZ(j),1)),1))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
        rectangle('Position',[2 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
    title(append('Neuron ',num2str(CC(ZZ(j),1))))
    if j==1 
        legend('Left','Right','Location','Northwest') 
    end

    end
    
end

end
end