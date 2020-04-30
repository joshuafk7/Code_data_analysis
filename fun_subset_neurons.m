function subset_neurons = fun_subset_neurons(neuron, subset)
    A = 1:size(neuron.C,1);
    B = setdiff(A,subset);
    neuron_temp = neuron;
    
    neuron_temp.A(:, B)=[];
    neuron_temp.C(B,:) =[];
    neuron_temp.C_raw(B,:)=[];
    neuron_temp.S(B,:) =[];
    neuron_temp.ids(B) =[];
    neuron_temp.tags(B) =[];
    subset_neurons = neuron_temp;
end