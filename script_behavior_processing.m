filename = dir('*.rhd');
taste_inputs = {'Suc','S_75','S_55','S_45','S_25','NaCl'};
direction_inputs = [1 1 1 2 2 2];
[tastes,unit, data,trial, summary] = process_intan_v4_behavior_only(filename, taste_inputs, direction_inputs);
DLC_output = dir('*.csv');
summary_w_DLC= fun_process_DLC_data(2, DLC_output, summary);