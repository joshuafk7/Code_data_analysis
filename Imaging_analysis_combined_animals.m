baseline =2;
comparison = 3;
taste_resp_neurons = fun_taste_responses_v2( tastes, combined_binning.binnedC6, combined_binning.binnedC7(:,:,2,:), baseline, comparison);
% plot_sig_neurons(taste_resp_neurons, plotbins, plotbinnedC2,summary)

comparison = 4;
delay_resp_neurons = fun_taste_responses_v2( tastes, combined_binning.binnedC6, combined_binning.binnedC7(:,:,2,:), baseline, comparison);
% plot_sig_neurons(taste_resp_neurons, plotbins, plotbinnedC2,summary)

comparison = 5;
choice_resp_neurons = fun_taste_responses_v2( tastes, combined_binning.binnedC6, combined_binning.binnedC7(:,:,2,:), baseline, comparison);
% plot_sig_neurons(choice_resp_neurons, plotbins, plotbinnedC2,summary)
%%
figure(726)

subplot(3,1,1)
popavgtaste1 = mean(zscore(squeeze(nanmean(combined_binning.plotbinnedC7(:,:,2,taste_resp_neurons),1)),0,1),2)';
hold on
plot(plotbins(1:end-1),popavgtaste1);
y=ylim;
rectangle('Position',[0 y(1) 2 y(2)-y(1)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[y(1) y(2)],'Color','r','LineStyle','--')
title('Taste Responsive Neurons')
% xlabel('Time (sec)')
ylabel('Spiking Z score')

subplot(3,1,2)
y=[];
popavgtaste2 = mean(zscore(squeeze(nanmean(combined_binning.plotbinnedC7(:,:,2,delay_resp_neurons),1)),0,1),2)';
hold on
plot(plotbins(1:end-1),popavgtaste2);
y=ylim;
rectangle('Position',[2 y(1) 2 y(2)-y(1)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[y(1) y(2)],'Color','r','LineStyle','--')
title('Delay Responsive Neurons')
% xlabel('Time (sec)')
ylabel('Spiking Z score')

subplot(3,1,3)
y=[];
popavgtaste3 = mean(zscore(squeeze(nanmean(combined_binning.plotbinnedC7(:,:,2,choice_resp_neurons),1)),0,1),2)';
hold on
plot(plotbins(1:end-1),popavgtaste3);
y=ylim;
rectangle('Position',[4 y(1) 2 y(2)-y(1)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[y(1) y(2)],'Color','r','LineStyle','--')
title('Choice Responsive Neurons')
% xlabel('Time (sec)')
ylabel('Spiking Z score')

%%
 p=[];
  t = categorical({'Suc','S\_75\_25','S\_55\_45','S\_45\_55','S\_25\_75','NaCl'});
 x=[100 75 55 45 25 0;];

%   t = categorical({'S\_55\_45','S\_45\_55'});
%  x=[ 55 45;];
  for i=1:length(taste_resp_neurons)
      taste_comparison=[];
  taste_comparison = squeeze(combined_binning.binnedC6(:,3,:,taste_resp_neurons(i)));
  [p,~,~] = anova1(taste_comparison,x,'off'); %compare all tastes in bin 3
  PP(i)=p;
  end
  idx1=find(PP<0.05);
  taste_selective=taste_resp_neurons(idx1);

%% clustering for neurons from 2way anova comparison
Q = squeeze(nanmean(combined_binning.binnedC6(:,3,:,taste_selective),1));
Q=zscore(Q);
QQ = linkage(Q','ward');
cutoff = median([QQ(end-1,3) QQ(end-0,3)]);
QQQ = cluster(QQ, 'Cutoff', cutoff, 'Depth',5,'Criterion','distance');
% QQQ = cluster(QQ, 'MaxClust', 8);
numclusters = length(unique(QQQ));
clust_neurons = NaN(length(taste_selective),numclusters);
for i=1:numclusters
    temp1 = taste_selective(find(QQQ == i));
    temp2=length(taste_selective)-length(temp1);
    temp3 = [temp1; NaN(temp2,1)];
   clust_neurons(:,i) =  temp3;
end

% % neuron_number = 1:length(neurons);
% clust1_neurons = taste_resp_2wayANOVA(find(QQQ == 1));
% clust2_neurons = taste_resp_2wayANOVA(find(QQQ == 2));
% clust3_neurons = taste_resp_2wayANOVA(find(QQQ == 3));
% clust4_neurons = taste_resp_2wayANOVA(find(QQQ == 4));
figure(716)
D=pdist(Q');
leaforder=optimalleaforder(QQ,D);
temp1=taste_selective;
for i=1:length(temp1)
    clustorder(i) = temp1(leaforder(i));
end
colors=rand(numclusters,3);
subplot(10,1,[1 2 3 4])
h=dendrogram(QQ,0,'ColorThreshold',cutoff,'Labels',num2str(taste_selective),'Reorder',leaforder);
% for i=1:length(h)
%     h(i).Color = colors(clustupdate(i),:);
%    
% end
% h(2).Color
set(gca,'fontsize',6)
title('Hierarchical Clustering - Taste Selective Neurons')
clustplot1 = zscore(squeeze(nanmean(combined_binning.binnedC6(:,3,:,clustorder,1))));
p=5;
j=2;
i=1;
c=categorical(clustorder, clustorder);
% for j=1:6
%     subplot(10,1,p)
%     hold on

p=5;


for j=1:6
subplot(10,1,p)
    ylabel(t(j))
    
    set(get(gca,'YLabel'),'Rotation',0)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    if j==6
       xlabel('Neurons') 
    end
for i=1:size(clustorder,2)
    hold on
%    if intersect(clustorder(i), clust_neurons(:,i)) ~=0
     
%       hold on
[~,y1] = find(clustorder(i)==clust_neurons);
    clustupdate(i)=y1;
    bar1= bar(c(i),clustplot1(j,i),'FaceColor',[colors(clustupdate(i),:)]);


end
p=p+1;
end



set(gcf, 'Renderer', 'painters');
%%
clusters = unique(clustupdate,'stable');
clustupdate=clustupdate';
% updated_clusters=NaN(length(taste_selective),numclusters);


for i=1:numclusters
     temp1=clustorder(find(clustupdate==clusters(i)));
    updated_clusters{1,i} = temp1;
end
for i=1:numclusters
    temp1 = updated_clusters{1,i};
    clustmean(:,i)=mean(zscore(squeeze(nanmean(combined_binning.binnedC6(:,3,:,temp1,1)))),2);
    temp1=[];
end

x=[100 75 55 45 25 0];
figure(667)
for i=1:numclusters
    subplot(2,1,i)
    scatter(x(1,:),clustmean(:,i)',40,colors(clusters(i),:),'filled')
    title(append('Cluster',num2str(i)))
    xlabel('Sucrose Concentration')
    ylabel('Spiking Z score')
end
set(gcf, 'Renderer', 'painters');
%%

pval=[];
temp=[];
CC=[];
GC=[];
q=1;
% GC=nchoosek(x,2);
comparison=3;
[pval,tastepairs] = fun_compare_tastes_v2(4,comparison, combined_binning.binnedC6,  tastes);
GC=tastepairs;
[CC(:,1),CC(:,2)]=find(pval<0.05);
for i=1:size(CC,1)
    if ~isempty(intersect(CC(i,1),taste_resp_neurons))
        temp(q,1) = CC(i,1);
        temp(q,2) = CC(i,2);
        q=q+1;
    end
end
CC = temp;
selective_neurons = unique(CC(:,1));
numpairs = [1:size(tastepairs,1)];
[GC(:,4),GC(:,3)] = groupcounts(CC(:,2)) ;

paircomp=nan(6,6);
for i=1:length(GC)
paircomp(GC(i,1),GC(i,2)) = GC(i,4);
end
figure(8312)
h6=heatmap(t,t,paircomp');
h6.MissingDataColor = [1 1 1];
h6.MissingDataLabel='';
h6.GridVisible='off';
h6.ColorbarVisible='off';
 set(gcf, 'Renderer', 'painters');
 %% linear fits
 
 comparison=3;
% x=[100; 55; 45; 0];
rsqu=[];
x=[100; 75; 55; 45; 25; 0];

for i=1:length(taste_resp_neurons)
y=[];
md1=[];
y = squeeze(nanmean(combined_binning.binnedC6(:,comparison,:,taste_resp_neurons(i)))); %find mean binnedC for each bin
% y=mean(zscore(squeeze(nanmean(combined_binning.binnedC6(:,3,:,taste_resp_neurons(i),1)))),2);
mdl = fitlm(x,y);
rsqu(i)=mdl.Rsquared.Adjusted;

end

mixture_coding=[];
mixture_coding = find(rsqu>.7 | rsqu<-.7);
length(mixture_coding)
length(find(rsqu>.7 ))
length(find(rsqu<-.7))
%%
figure(736)

for i=1:length(mixture_coding)
    subplot(4,4,i)
    y = squeeze(nanmean(combined_binning.binnedC(:,comparison,:,taste_resp_neurons(mixture_coding(i)))));
    p = polyfit(x,y,1); 
    f = polyval(p,x); 
    plot(x,y,'o',x,f,'-') 
    title(append('Neuron ',num2str(taste_resp_neurons(mixture_coding(i)))))
% legend('data','linear fit')
end