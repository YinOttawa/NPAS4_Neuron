%% to calculate the peak during RB
%to calculate the eligible peak
temp_dff = [];
peak_ind = struct;
peak_ind_2ndmean = struct;

for n = 1:length(df_Ff_success)
    temp_dff = df_Ff_success(n).reaching_dff;
    [mm nn] = size(temp_dff);
    dff_cell = [];    
    
    NP = 1;
    NN = 1;
    for z = 1:mm
        dff_cell = temp_dff(z,:);
        threshold_cell = df_baseline_threshold(z);
        dff_cell_bina = dff_cell;
        dff_cell_bina(dff_cell >= threshold_cell) = 1;
        dff_cell_bina(dff_cell < threshold_cell) = 0;

        peak_count = 1;
        peak_start = 0;
        Peak_end =0;
             
        peak_search = 1;
        while peak_search <= length(dff_cell_bina)-3
            bina_rest = [];
            if dff_cell_bina(peak_search) == 1
                peak_search = peak_search+1;
                if dff_cell_bina(peak_search) == 1
                    peak_search = peak_search+1;
                    if dff_cell_bina(peak_search) == 1                        
                        peak_start = peak_search-2;
                        
                        % find three 0 in a row, defined by baseline_mean
                        % instead of baseline_threshold
                        dff_cell = temp_dff(z,:);
                        threshold_cell_2nd = df_baseline_mean(z);
                        dff_cell_bina_mean = dff_cell;
                        dff_cell_bina_mean(dff_cell >= threshold_cell_2nd) = 1;
                        dff_cell_bina_mean(dff_cell < threshold_cell_2nd) = 0;

                        bina_rest = dff_cell_bina_mean(peak_search+1:end);  
                        find_zero = find(~bina_rest);% locate the  0
                        find_zero_diff = diff(find_zero);                        
                        find_zero_diff(find_zero_diff>1) = 0;
                        zero_count = 0;
                        for count_n = 1:length(find_zero_diff)
                            if find_zero_diff(count_n) ==0
                                zero_count = 0;
                            else
                                zero_count = zero_count+1;
                                if zero_count ==2  % here are three 0 in a row
                                    temp_end = find_zero(count_n-1); % locate the first zero (of three 0 in a row)
                                    peak_end = peak_search + temp_end-1; %
                                    break
                                end
                            end
                        end
                        
                        % if zero_count <3, then the whole are an event
                        if zero_count <2
                            peak_end = length(dff_cell_bina);
                        end                  
                        
                        
                        peak_ind_2ndmean(n).trial_num = n;
                        peak_ind_2ndmean(n).frames = nn;                        
                        peak_ind_2ndmean(n).peak(z).cell_num = z;
                        peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_count = peak_count;
                        peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_start = peak_start;
                        peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_end = peak_end;
                        peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_duration = peak_end - peak_start +1;
                      

                        dff_event = dff_cell(peak_start:peak_end);
                        peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_max_dff = max(dff_event);
                        
                        peak_count = peak_count +1;
                        peak_search = peak_end+1;
%                         
                    else
                        peak_search = peak_search+1;
                    end                    
                else
                    peak_search = peak_search+1;
                end
                
            else
                peak_search = peak_search+1;
            end
            

            
            
            if peak_start == 0;
                peak_ind_2ndmean(n).trial_num = n;
                peak_ind_2ndmean(n).frames = nn;
                peak_ind_2ndmean(n).peak(z).cell_num = z;                
                peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_count = NaN;
                peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_start = NaN;
                peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_end = NaN;
                peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_duration = NaN;
                peak_ind_2ndmean(n).peak(z).peak(peak_count).peak_max_dff = NaN;
                

            end
            
            % to calculate the width and frequency
            cal_table = struct2table(peak_ind_2ndmean(n).peak(z).peak);
            cal_mean = nanmean(cal_table.peak_duration);
            cal_max = nanmax(cal_table.peak_duration); 
            cal_dff_max = nanmax(cal_table.peak_max_dff);
            cal_dff_max_mean = nanmean(cal_table.peak_max_dff);
            peak_ind_2ndmean(n).peak(z).duration = cal_mean;  
            peak_ind_2ndmean(n).peak(z).duration_max = cal_max;  
            peak_ind_2ndmean(n).peak(z).duration_sec = cal_mean/frame_rate;
            peak_ind_2ndmean(n).peak(z).duration_max_sec = cal_max/frame_rate;
            peak_ind_2ndmean(n).peak(z).peak_max_dff = cal_dff_max;
            peak_ind_2ndmean(n).peak(z).peak_max_dff_mean = cal_dff_max_mean;
            
            if isnan(cal_mean)
                peak_ind_2ndmean(n).peak(z).events = NaN;
                peak_ind_2ndmean(n).peak(z).frequency_sec = NaN;
            else               
                peak_ind_2ndmean(n).peak(z).events = length(peak_ind_2ndmean(n).peak(z).peak);
                peak_ind_2ndmean(n).peak(z).frequency_sec = peak_ind_2ndmean(n).peak(z).events/nn*frame_rate;
            end
            
            if ismember(z, NRAM_pos_ind)
                peak_ind_2ndmean(n).peak(z).NRAM_pos = 1;                
            else
                peak_ind_2ndmean(n).peak(z).NRAM_pos = 0;
            end

        end      
    end
end


% plots
duration_total = [];
duration_max_total = [];
frequency_total = [];
peak_max_dff_total = [];
peak_max_dff_mean_total = [];
for xx = 1: length(peak_ind_2ndmean)
    cal_table = [];
    cal_table = struct2table(peak_ind_2ndmean(xx).peak);    
    duration_total = [duration_total cal_table.duration_sec];
    duration_max_total = [duration_max_total cal_table.duration_max_sec];
    frequency_total = [frequency_total cal_table.frequency_sec];
    peak_max_dff_total = [peak_max_dff_total cal_table.peak_max_dff];
    peak_max_dff_mean_total = [peak_max_dff_mean_total cal_table.peak_max_dff_mean];
end

duration_mean_pos = nanmean(duration_total(NRAM_pos_ind,:))';
duration_mean_neg = nanmean(duration_total(NRAM_neg_ind,:))';

duration_max_pos = nanmean(duration_max_total(NRAM_pos_ind,:))';
duration_max_neg = nanmean(duration_max_total(NRAM_neg_ind,:))';
      
frequency_mean_pos = nanmean(frequency_total(NRAM_pos_ind,:))';
frequency_mean_neg = nanmean(frequency_total(NRAM_neg_ind,:))';

peak_max_dff_pos = nanmean(peak_max_dff_total(NRAM_pos_ind,:))';
peak_max_dff_neg = nanmean(peak_max_dff_total(NRAM_neg_ind,:))';

peak_max_dff_mean_pos = nanmean(peak_max_dff_mean_total(NRAM_pos_ind,:))';
peak_max_dff_mean_neg = nanmean(peak_max_dff_mean_total(NRAM_neg_ind,:))';

data_RB_calcium(session_num).event_2ndmean.peak_ind = peak_ind_2ndmean;
data_RB_calcium(session_num).event_2ndmean.duration_mean_pos = duration_mean_pos ;
data_RB_calcium(session_num).event_2ndmean.duration_mean_neg =duration_mean_neg;
data_RB_calcium(session_num).event_2ndmean.duration_max_pos = duration_max_pos;
data_RB_calcium(session_num).event_2ndmean.duration_max_neg = duration_max_neg;      
data_RB_calcium(session_num).event_2ndmean.frequency_mean_pos = frequency_mean_pos;
data_RB_calcium(session_num).event_2ndmean.frequency_mean_neg = frequency_mean_neg;
data_RB_calcium(session_num).event_2ndmean.peak_max_dff_pos = peak_max_dff_pos;
data_RB_calcium(session_num).event_2ndmean.peak_max_dff_neg = peak_max_dff_neg;
data_RB_calcium(session_num).event_2ndmean.peak_max_dff_pos = peak_max_dff_mean_pos;
data_RB_calcium(session_num).event_2ndmean.peak_max_dff_neg = peak_max_dff_mean_neg;

        
% plot for probability
    prob_duration_pos = [];
    prob_duration_neg = [];
    prob_max_dff_pos = [];
    prob_max_dff_neg = [];
    prob_freq_pos = [];
    prob_freq_neg = [];
for n = 1:length(df_Ff_success)  % trial number
    temp_dff = df_Ff_success(n).reaching_dff;
    [mm nn] = size(temp_dff);
%     dff_cell = [];     

    for z = 1:mm   % cell number  
        if ismember(z,NRAM_pos_ind)            
            cal_table2 = [];
            cal_table2 = struct2table(peak_ind_2ndmean(n).peak(z).peak);        
            prob_duration_pos = [prob_duration_pos;cal_table2.peak_duration];
            prob_max_dff_pos = [prob_max_dff_pos;cal_table2.peak_max_dff];
            prob_freq_pos = [prob_freq_pos;peak_ind_2ndmean(n).peak(z).frequency_sec];
        else
            cal_table2 = [];
            cal_table2 = struct2table(peak_ind_2ndmean(n).peak(z).peak);        
            prob_duration_neg = [prob_duration_neg;cal_table2.peak_duration];
            prob_max_dff_neg = [prob_max_dff_neg; cal_table2.peak_max_dff];
            prob_freq_neg = [prob_freq_neg; peak_ind_2ndmean(n).peak(z).frequency_sec];
        end
    end
end

prob_duration_pos = prob_duration_pos(~isnan(prob_duration_pos)); %% delete NaN
prob_duration_neg = prob_duration_neg(~isnan(prob_duration_neg));
prob_max_dff_pos = prob_max_dff_pos(~isnan(prob_max_dff_pos));
prob_max_dff_neg = prob_max_dff_neg(~isnan(prob_max_dff_neg));
prob_freq_pos = prob_freq_pos(~isnan(prob_freq_pos));
prob_freq_neg = prob_freq_neg(~isnan(prob_freq_neg));


%plot duration
figure
x_pos_prob = sort(prob_duration_pos)/frame_rate;
y_pos_prob = [1:length(prob_duration_pos)]/length(prob_duration_pos);
plot(x_pos_prob, y_pos_prob,'r')

hold on
x_neg_prob = sort(prob_duration_neg)/frame_rate;
y_neg_prob = [1:length(prob_duration_neg)]/length(prob_duration_neg);
plot(x_neg_prob, y_neg_prob,'k')
title_text = ['2nd mean event Duration'];
title(title_text);
hold off

%plot max_dff
figure
x_pos_prob = sort(prob_max_dff_pos);
y_pos_prob = [1:length(prob_max_dff_pos)]/length(prob_max_dff_pos);
plot(x_pos_prob, y_pos_prob,'r')

hold on
x_neg_prob = sort(prob_max_dff_neg);
y_neg_prob = [1:length(prob_max_dff_neg)]/length(prob_max_dff_neg);
plot(x_neg_prob, y_neg_prob,'k')
title_text = ['2nd mean event max dff'];
title(title_text);
        
% plot freq        
figure
x_pos_prob = sort(prob_freq_pos);
y_pos_prob = [1:length(prob_freq_pos)]/length(prob_freq_pos);
plot(x_pos_prob, y_pos_prob,'r')

hold on
x_neg_prob = sort(prob_freq_neg);
y_neg_prob = [1:length(prob_freq_neg)]/length(prob_freq_neg);
plot(x_neg_prob, y_neg_prob,'k')
title_text = ['2nd mean event frequency'];
title(title_text);
                


data_RB_output = struct;
data_RB_output(1).ITI_RB_dff = nanmean(data_RB_calcium(session_num).ITI_dff_mean.ITI_dff_mean_total);
data_RB_output(2).ITI_RB_dff = nanmean(data_RB_calcium(session_num).reaching_dff_mean_active.reaching_dff_mean_total_a);

data_RB_output(1).ITI = nanmean(data_RB_calcium(session_num).ITI_dff_mean.ITI_dff_mean_neg);
data_RB_output(2).ITI = nanmean(data_RB_calcium(session_num).ITI_dff_mean.ITI_dff_mean_pos);

data_RB_output(1).mean_dff_pos_vs_neg = nanmean(data_RB_calcium(session_num).reaching_dff_mean_active.reaching_dff_mean_neg_a);
data_RB_output(2).mean_dff_pos_vs_neg = nanmean(data_RB_calcium(session_num).reaching_dff_mean_active.reaching_dff_mean_pos_a);

data_RB_output(1).Prop_active = nanmean(data_RB_calcium(session_num).prop_active.prop_active_neg);
data_RB_output(2).Prop_active = nanmean(data_RB_calcium(session_num).prop_active.prop_active_pos);

data_RB_output(1).reliability = nanmedian(data_RB_calcium(session_num).reliability.prop_rel_neg);
data_RB_output(2).reliability = nanmedian(data_RB_calcium(session_num).reliability.prop_rel_pos);

data_RB_output(1).event_duration = nanmean(data_RB_calcium(session_num).event_2ndmean.duration_mean_neg);
data_RB_output(2).event_duration = nanmean(data_RB_calcium(session_num).event_2ndmean.duration_mean_pos);

data_RB_output(1).event_max_dff = nanmean(data_RB_calcium(session_num).event_2ndmean.peak_max_dff_neg);
data_RB_output(2).event_max_dff = nanmean(data_RB_calcium(session_num).event_2ndmean.peak_max_dff_pos);


dff_mean_cell = [];
for zzzz = 1:size(df_Ff_success,2)
    dff_mean_cell = [dff_mean_cell df_Ff_success(zzzz).reaching_dff_mean];
end

data_RB_calcium(session_num).reaching_dff_mean_active.reaching_dff_mean_neg_cell = dff_mean_cell(NRAM_neg_ind_a);
data_RB_calcium(session_num).reaching_dff_mean_active.reaching_dff_mean_pos_cell = dff_mean_cell(NRAM_pos_ind_a);


data_RB_output(1).mean_dff_pos_vs_neg_cell = nanmedian(dff_mean_cell(NRAM_neg_ind_a));
data_RB_output(2).mean_dff_pos_vs_neg_cell = nanmedian(dff_mean_cell(NRAM_pos_ind_a));

zzzz_dff_mean_cell = mean(dff_mean_cell,2);
xxxx_new_order = sortrows(cell_map_all,2); %sort based on the 2nd column

    
