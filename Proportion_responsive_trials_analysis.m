numTrials = length(trial);
normResponses = zeros(size(binnedC2, 3),size(binnedC2, 2));
numBins = 1:(length(bins));
%find percentage of trials with responses in a given bin for all neurons
for i =1:size(binnedC2, 3)
    col=[];
    [~,col]=find(binnedC2(:,:,i)); %column index corresponds to the number of non zero bins over all trials
    normResponses(i,:) = histcounts(col,numBins)./numTrials;% fraction of trials with responses in each bin
end
shuffled=[];
for j = 1:1000
    shuffled = zeros(size(binnedC2));
    for i =1:length(trial)
        randBins = randperm(20);
        shuffled(i,:,:) = binnedC2(i,randBins, :);
    end
    for i =1:size(shuffled, 3)
        col=[];
        [~,col]=find(shuffled(:,:,i)); %column index corresponds to the number of non zero bins over all trials
        normResponsesShuff(i,:,j) = histcounts(col,numBins)./numTrials;% fraction of trials with responses in each bin
    end
end
for i = 1:size(normResponsesShuff, 1)

temp=[];
temp = squeeze(normResponsesShuff(i, :,:))';
% sem=std(temp)/sqrt(1000);
twostdevs(i,:)=std(temp).*3;

end
p=1;
modulatedneurons = [];
for i = 1:size(normResponsesShuff, 1)
    if~isempty(find(normResponses(i,8:end)>twostdevs(i,8:end)))
       modulatedneurons(p) = i;
       p=p+1;
    end
end

neuron2pick = 6;
temp=[];
temp = squeeze(normResponsesShuff(neuron2pick, :,:))';
figure
subplot(2,1,1)
bar(bins(1:end-1),normResponses(neuron2pick,:))
subplot(2,1,2)
bar(bins(1:end-1),mean(temp))
hold on
errorbar(bins(1:end-1), mean(temp),twostdevs(neuron2pick,:)) 


[M,I] = max(normResponses');
length(histcounts(col,numBins))