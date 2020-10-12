
size(binnedC_JK146)
combined_trials = size(binnedC_JK141,[1 ])+size(binnedC_JK146,[1 ]);
bin_num = size(binnedC_JK141,2);
taste_num = size(binnedC_JK141,3);
combined_neurons = size(binnedC_JK141,[4 ])+size(binnedC_JK146,[4 ]);

binnedC_combined = NaN(combined_trials,bin_num,taste_num,combined_neurons);

binnedC_combined(1:size(binnedC_JK141,1),:,:,1:size(binnedC_JK141,4)) = binnedC_JK141;
binnedC_combined((size(binnedC_JK141,1)+1):(size(binnedC_JK141,1))+size(binnedC_JK146,1),:,:,(size(binnedC_JK141,4)+1):(size(binnedC_JK141,4))+size(binnedC_JK146,4)) = binnedC_JK146;
size(binnedC_combined)

binned_names = fieldnames(JK141);

for i = 1:length(binned_names)
   if contains( binned_names(i),'binned')
       combined_binning.(binned_names{i}) = [];
      
   end
    
end


binned_names_only = fieldnames(combined_binning);


JK141_size = size(JK141.(binned_names_only{1}));
JK145_size = size(JK145.(binned_names_only{1}));
JK146_size = size(JK146.(binned_names_only{1}));



for i = 1:length(binned_names_only)
    individual_size = [size(JK141.(binned_names_only{i}));size(JK145.(binned_names_only{i}));size(JK146.(binned_names_only{i}))];
%    combined_size = zeros(1,size(individual_size,2));
   combined_size = individual_size(1,:);
    for j = 1:length(combined_size)
      
        if j == 1
            combined_size(j) = sum(individual_size(:,j));
        end
        if j == length(combined_size)
           combined_size(j) = sum(individual_size(:,j)) ;
        end
       
    end
    combined_binning.(binned_names_only{i}) = NaN(combined_size);
    cumulative_size = cumsum(individual_size);
    for p = 1:size(individual_size,1)
    for j = 1:size(individual_size,2)
        if j ==1 || j==size(individual_size,2)
        if p==1
        dims{p,j} = 1:individual_size(p,j);
        end
        if p>1
        dims{p,j} = cumulative_size(p-1,j)+1:cumulative_size(p,j);
        end
        else
          dims{p,j} = 1:individual_size(p,j);  
        end
    end
    end
    
    combined_binning.(binned_names_only{i})(dims{1,:}) = JK141.(binned_names_only{i});
    combined_binning.(binned_names_only{i})(dims{2,:}) = JK145.(binned_names_only{i});
    combined_binning.(binned_names_only{i})(dims{3,:}) = JK146.(binned_names_only{i});
    individual_size=[];
    combined_size=[];
    dims=[];
end