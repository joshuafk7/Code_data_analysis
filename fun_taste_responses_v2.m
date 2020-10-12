function sig_neurons = fun_taste_responses_v2( tastes, binnedC, binnedC2, baseline, comparison)
pval3=NaN(length(tastes),size(binnedC,4));
x=[];
idx=1;
for j =1:length(tastes) %checking for each taste individually
    for i=1:size(binnedC,4)
        %pval - dim1 = taste, dim2 = neuron, dim3 = bin under comparison to
        %baseline
        pval3(j,i)=ranksum(binnedC(:,baseline,j,i),binnedC(:,comparison,j,i));
    end
end
x=[];
[x(:,1),x(:,2)]=find(pval3<(0.05)); %find neurons where pvalue is less than 0.05
sig_neurons=[];
sig_neurons = unique(x(:,2)); %List of neurons with significant p value (some may respond to multiple tastes)
temp1 = squeeze(nanmean(binnedC2(:,baseline,sig_neurons)));
temp2= squeeze(nanmean(binnedC2(:,comparison,sig_neurons)));
larger = temp1<temp2;
sig_neurons = sig_neurons(larger); %only choose excitatory responses
fprintf(append('Number of responsive neurons: ', num2str(length(sig_neurons)), '/',num2str(size(binnedC,4)),'\n'))

end