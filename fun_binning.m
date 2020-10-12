%% this function creates a  binned matrix for all neurons/conditions, 
%if n = 1, trial x bin x tastes x neuron
%if n=2, trial x bin x neuron
%n=3 trial x bin x correct/error x neuron
%since only one taste is delivered per trial, all other trials are NaN
%

function binnedS = fun_binning(n,trial, bins, tastes, neurons)
switch n
    case 1
        binnedS = NaN(length(trial),length(bins)-1,length(tastes),length(neurons));
        for z = 1:length(neurons)
            x=1;
            for p=1:length(tastes)
                for i=1:length(neurons(z).(tastes{p}).spikes)
                    for j=1:length(bins)-1
                        idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                        binnedS(x,j,p,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx));
                    end
                    x=x+1;
                end
            end
        end
    case 2
        binnedS = NaN(length(trial),length(bins)-1,length(neurons));
        for z = 1:length(neurons)
            x=1;
            for p=1:length(tastes)
                for i=1:length(neurons(z).(tastes{p}).spikes)
                    for j=1:length(bins)-1
                        idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                        binnedS(x,j,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx));
                    end
                    x=x+1;
                end
            end
        end
        case 3
        binnedS = NaN(length(trial),length(bins)-1,2,length(neurons));
        for z = 1:length(neurons)
            x=1;
            for p=1:length(tastes)
                for i=1:length(neurons(z).(tastes{p}).spikes)
                    if cell2mat(neurons(z).(tastes{p}).correct_choice(i,1)) == 1
                    r=cell2mat(neurons(z).(tastes{p}).L_R_trial);
                    else
                        r = 2+ cell2mat(neurons(z).(tastes{p}).L_R_trial);
                    end
                    for j=1:length(bins)-1
                        idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                        binnedS(x,j,r,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx));
                    end
                    x=x+1;
                end
            end
        end
        case 4
        binnedS = NaN(length(trial),length(bins)-1,length(tastes)*2,length(neurons));
        for z = 1:length(neurons)
            x=1;
            for p=1:length(tastes)
                for i=1:length(neurons(z).(tastes{p}).spikes)
                    if cell2mat(neurons(z).(tastes{p}).correct_choice(i,1)) == 1
                    r=p;
                    else
                        r = length(tastes)+p;
                    end
                    for j=1:length(bins)-1
                        idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                        binnedS(x,j,r,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx));
                    end
                    x=x+1;
                end
            end
        end
    case 5
     for z = 1:length(neurons)
            x=1;
            for p=1:length(tastes)
                for i=1:length(neurons(z).(tastes{p}).spikes)
                    
                    r=cell2mat(neurons(z).(tastes{p}).L_R_trial);
                    
                    
                    for j=1:length(bins)-1
                        idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                        binnedS(x,j,r,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx));
                    end
                    x=x+1;
                end
            end
     end
        case 7
        binnedS = NaN(length(trial),length(bins)-1,2,length(neurons));
        for z = 1:length(neurons)
            x=1;
            for p=1:length(tastes)
                for i=1:length(neurons(z).(tastes{p}).spikes)
                    r=cell2mat(neurons(z).(tastes{p}).correct_choice(i,1))+1;
                   
                    for j=1:length(bins)-1
                        idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                        binnedS(x,j,r,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx));
                    end
                    x=x+1;
                end
            end
        end
end
end