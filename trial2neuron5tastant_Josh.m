%% this code is to reorganize the trial structure into neuron structure

function neuron = trial2neuron5tastant(trial,tastes);
% tastes = {'T_1','T_2','T_3','T_4','T_5','T_6','T_7'};
for n = 1:size(trial(1).traces,1)
    h = ones(1,length(tastes)); %each position in h keeps track of number of trials for each taste
    for i = 1: length(trial)
        for q =1:length(tastes)
            if ~isnan(trial(i).(tastes{q}))
                
                
                neuron(n).(tastes{q}).spikes{h(q),1}= trial(i).spikes(n,:); %grab spike values for each trial for that neuron
                neuron(n).(tastes{q}).traces{h(q),1}= trial(i).traces(n,:); %grab trace values for each trial for that neuron
                neuron(n).(tastes{q}).Frames{h(q),1}=trial(i).Frames_adjusted; %grab frames aka timestamps for each frame above
                neuron(n).(tastes{q}).centSp{h(q),1} = trial(i).centSp; %cent sp licks
                neuron(n).(tastes{q}).LeftSp{h(q),1} = trial(i).LeftSp; %right sp licks
                neuron(n).(tastes{q}).RightSp{h(q),1} = trial(i).RightSp; %left sp licks
                neuron(n).(tastes{q}).L_R_trial{h(q),1} = trial(i).L_R_trial; %left or right trial
                neuron(n).(tastes{q}).correct_choice{h(q),1} = trial(i).correct_choice; %correct or error
                neuron(n).(tastes{q}).total_performance{h(q),1} = trial(i).total_performance; %total cumulative performance until that trial
                neuron(n).(tastes{q}).(append((tastes{q}),'_performance')){h(q),1} = trial(i).(append((tastes{q}),'_performance'));
                %                 T_1_licks{a,1}    = trial(i).centSp;
                %                 T_1_licks{a,2}    = trial(i).RightSp;
                %                 T_1_licks{a,3}    = trial(i).LeftSp;
                h(q) = h(q)+1;
            end
        end
    end
end

neuron = orderfields(neuron,tastes);

end


    





%%

%
%% statistical test

