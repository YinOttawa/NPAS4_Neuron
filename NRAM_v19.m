close all
clear
clc

% mouse name
mouse_name = 'XY312';
% training day
T_Day = 'PR10';
csv_path = 'C:\Users\xyin2\Desktop\XY_Science_2021\csv_data_DLC\';
plane1_dff_path = 'F:\Suite2P_output\dff\';
aligned_dff_path = 'F:\Suite2P_analysis\';  %df_Ff_new is the aligned dff data (calcium data before first tone were removed)
dff_file = [plane1_dff_path mouse_name '_plane1_dff.mat'];
load(dff_file)

frame_rate = 30.03;
threshold_per = 0.5;  % 50% threshold
% %%%%%%%%%%%%%
% ITI_dura = 10;  % changable
Fail_dura = 6; % changable
% %%%%%%%%%%%%%
% ITI_frame = round(ITI_dura*frame_rate);

ITI_start = 11; %11 seconds before tone
ITI_end = 5 ;% 6 seconds duration (=11-5)
ITI_frame = round(ITI_start*frame_rate); %ITI start
ITI_frame_end = round(ITI_end*frame_rate); %ITI start

%% calculate baseline and active threshold
% find the session in the matrix, since PR01 in line1; PR04 in line2 etc.
session_all = {'PR01' 'PR04'  'PR07'  'PR10'};
session_day = {T_Day};
x_ind = strcmp(session_all,session_day);  %% compare strings
session_num = find(x_ind, 1);

df = data(session_num).df_Ff;
[num_cell num_frame] = size(df);

num_spikes = [];
df_baseline_sd = 0;
df_baseline_mean_individual = 0;
df_baseline_threshold = 0;
df_sort = [];
y = [];
active_cell_session = zeros(1, num_cell);
active_cell_num = [];
ii = 1;
% calculate threshold by using 50% of sorted dff
for curr_cell = 1:num_cell    
    % sort df_curr_cell   
    df_curr = df(curr_cell,:);
    df_sort = sort(df_curr);
    df_baseline = df_sort(1:num_frame*threshold_per);    
    df_baseline_sd = std(df_baseline);
    df_baseline_mean_individual = mean(df_baseline,2);
    df_baseline_threshold(curr_cell,:) = df_baseline_mean_individual + 3*df_baseline_sd; 
    df_baseline_mean(curr_cell,:) = df_baseline_mean_individual;
 
    %find spike event (3 consecutive frames over threshold)
    bina_dff_sess = ones(1,num_frame);    
    for x = 1: num_frame
        if df_curr(x)< df_baseline_threshold(curr_cell)
            bina_dff_sess(x) = 0;
        end
    end
    frame_count = 0;
    for frame_c = 1: num_frame        
        if bina_dff_sess(frame_c) == 0
            frame_count = 0;
        else 
            frame_count = frame_count + 1;
            if frame_count == 3
                active_cell_session(curr_cell)=1;
                active_cell_num(ii) = curr_cell;  %% to indicate which cell is active
                ii = ii+1;
                break
            end
        end
    end       
end


clearvars -except df_baseline_threshold df_baseline_mean active_cell_session active_cell_num mouse_name T_Day csv_path dff_file frame_rate new_cell_index_map session_num data aligned_dff_path ITI_frame ITI_frame_end Fail_dura
%% calculate the active cells
dff_aligned_file = [aligned_dff_path mouse_name '\' mouse_name '_' T_Day '\' mouse_name '_' T_Day '_mean_df_Ff.mat'];
load(dff_aligned_file)

for num_trial = 1:size(df_Ff_success, 2)
    % locate the reaching bout dff
    start_char = data_video_success(num_trial).reaching_start_sec; 
    temp = strsplit(start_char,':');
    start_sec = str2double(temp(2));
    
    end_char = data_video_success(num_trial).reaching_end_sec;
    temp = strsplit(end_char,':');
    end_sec = str2double(temp(2));

    start_ind = round((start_sec + 1)*frame_rate);  %%+1 because the first second is tone time, the cutted video is from (tone_start+1)
    end_ind = round((end_sec + 1)*frame_rate); %%+1 because the first second is tone time
    reaching_dff = df_Ff_success(num_trial).df_Ff(:,start_ind:end_ind);
    df_Ff_success(num_trial).reaching_dff = reaching_dff;
    df_Ff_success(num_trial).reaching_dff_mean = mean(reaching_dff,2);
    
    % find spike event (3 consecutive frames over threshold)
    [n_cell m_frame] = size(reaching_dff);      
    bina_dff = [];
    temp1 = [];
    active_cell = zeros(1, n_cell);
    for num_cell = 1:n_cell
        % binarize the dff to make the analysis easier
        bina_dff = ones(1,m_frame);
        temp1 = reaching_dff(num_cell,:);
        for num_frame = 1: m_frame
            if temp1(num_frame)< df_baseline_threshold(num_cell)
                bina_dff(num_frame) = 0;
            end
        end
        frame_count = 0;
        for frame_c = 1: m_frame        
            if bina_dff(frame_c) == 0
                frame_count = 0;
            else 
                frame_count = frame_count + 1;
                if frame_count == 3
                    active_cell(num_cell)=1;
                    break
                end
            end
        end        
    end
    df_Ff_success(num_trial).active_cell = active_cell';         
end


data_RB_calcium = [];

%% read NRAM som cells and load the cell_map from session to session
file_NRAM = [csv_path mouse_name '\' mouse_name '_' T_Day '_NRAM_SOM.xlsx'];
file_NRAM_pos = readtable(file_NRAM);
NRAM_pos = table2array(file_NRAM_pos);

cell_map_all = [];
for mm = 1:session_num;  %%session4 = 2 (session_num was calculated earlier)
    cell_map_all = [cell_map_all; new_cell_index_map(mm).sort_cells];
end
cell_map_curr = cell_map_all(:,session_num);

cell_ind = [];
for num_som = 1:length(NRAM_pos)
    cell_ind(num_som) = cell_map_curr(NRAM_pos(num_som));  %%triple checked, it is right   
%     cell_ind(num_som) = find(cell_map_curr==NRAM_pos(num_som));
end
NRAM_pos_ind = sort(cell_ind);
for i = 1:length(NRAM_pos_ind)  %% exclude the non-active cells during the whole session
    curr_cell = NRAM_pos_ind(i);
    if ismember(curr_cell,active_cell_num) == 0
        NRAM_pos_ind(i)=NaN;
    end
end
% NRAM_pos_ind_low_index = NRAM_pos_ind;
NRAM_pos_ind = NRAM_pos_ind(~isnan(NRAM_pos_ind)); %% delete NaN
    
temp_cell_pool = [1:1:n_cell];
temp_cell_pool(NRAM_pos_ind) = []; %%remove the NRAM_pos_ind
NRAM_neg_ind = temp_cell_pool;
for i = 1:length(NRAM_neg_ind)  %% exclude the non-active cells during the whole session
    curr_cell = NRAM_neg_ind(i);
    if ismember(curr_cell,active_cell_num) == 0
        NRAM_neg_ind(i)=NaN;
    end
end
NRAM_neg_ind = NRAM_neg_ind(~isnan(NRAM_neg_ind)); %% delete NaN

%% to calculate the reliability
active_rel = [];
for num_trial = 1:size(df_Ff_success, 2)
    active_rel = [active_rel df_Ff_success(num_trial).active_cell];
end
prop_rel = sum(active_rel,2)/size(active_rel,2);
prop_rel_N_pos = prop_rel(NRAM_pos_ind);
prop_rel_N_neg = prop_rel(NRAM_neg_ind);

prop_rel_N_pos = sort(prop_rel_N_pos);
prop_rel_N_neg = sort(prop_rel_N_neg);

x_pos = [1:1:length(prop_rel_N_pos)]/length(prop_rel_N_pos);
x_neg = [1:1:length(prop_rel_N_neg)]/length(prop_rel_N_neg);

data_RB_calcium(session_num).reliability.prop_rel_pos = prop_rel_N_pos;
data_RB_calcium(session_num).reliability.prop_rel_neg = prop_rel_N_neg;
% figure
% plot(prop_rel_N_pos, x_pos);
% hold on
% plot(prop_rel_N_neg, x_neg );




%% to calculate active cell proportion from each trial
for num_trial = 1:size(df_Ff_success, 2)
    active_cell = [];
    active_cell = df_Ff_success(num_trial).active_cell;
    prop_active_pos(num_trial) = sum(active_cell(NRAM_pos_ind))/length(NRAM_pos_ind);
    prop_active_neg(num_trial) = sum(active_cell(NRAM_neg_ind))/length(NRAM_neg_ind); 
end

prop_active_pos_s = sort(prop_active_pos);
prop_active_neg_s = sort(prop_active_neg);

data_RB_calcium(session_num).prop_active.prop_active_pos = prop_active_pos';
data_RB_calcium(session_num).prop_active.prop_active_neg = prop_active_neg';

% figure
% plot(prop_active_pos_s);
% hold on
% plot(prop_active_neg_s);


%% to calculate mean during reaching bout
for num_trial = 1:size(df_Ff_success, 2)
    reaching_dff_mean = [];
    reaching_dff_mean = df_Ff_success(num_trial).reaching_dff_mean;
    reaching_dff_mean_pos(num_trial) = nanmean(reaching_dff_mean(NRAM_pos_ind));
    reaching_dff_mean_neg(num_trial) = nanmean(reaching_dff_mean(NRAM_neg_ind));
    reaching_dff_mean_total(num_trial) = nanmean(reaching_dff_mean); 
end
reaching_dff_mean_pos_s = sort(reaching_dff_mean_pos);
reaching_dff_mean_neg_s = sort(reaching_dff_mean_neg);

data_RB_calcium(session_num).reaching_dff_mean.reaching_dff_mean_pos = reaching_dff_mean_pos';
data_RB_calcium(session_num).reaching_dff_mean.reaching_dff_mean_neg = reaching_dff_mean_neg';
data_RB_calcium(session_num).reaching_dff_mean.reaching_dff_mean_total = reaching_dff_mean_total';

% figure
% plot(reaching_dff_mean_pos_s);
% hold on
% plot(reaching_dff_mean_neg_s);   

%% to calculate mean during reaching bout, only active cells
for num_trial = 1:size(df_Ff_success, 2)
    reaching_dff_mean = [];
    reaching_dff_mean = df_Ff_success(num_trial).reaching_dff_mean;  
    reaching_dff_mean_active = reaching_dff_mean;
    active_ind = df_Ff_success(num_trial).active_cell;
    for num_cell = 1:length(reaching_dff_mean)
        if active_ind(num_cell) == 0
            reaching_dff_mean_active(num_cell) = NaN;
        end
    end           
    reaching_dff_mean_pos_a(num_trial) = nanmean(reaching_dff_mean_active(NRAM_pos_ind));
    reaching_dff_mean_neg_a(num_trial) = nanmean(reaching_dff_mean_active(NRAM_neg_ind));
    reaching_dff_mean_total_a(num_trial) = nanmean(reaching_dff_mean_active);
end
reaching_dff_mean_pos_s_a = sort(reaching_dff_mean_pos_a);
reaching_dff_mean_neg_s_a = sort(reaching_dff_mean_neg_a);

data_RB_calcium(session_num).reaching_dff_mean_active.reaching_dff_mean_pos_a = reaching_dff_mean_pos_a';
data_RB_calcium(session_num).reaching_dff_mean_active.reaching_dff_mean_neg_a = reaching_dff_mean_neg_a';
data_RB_calcium(session_num).reaching_dff_mean_active.reaching_dff_mean_total_a = reaching_dff_mean_total_a';
% figure
% plot(reaching_dff_mean_pos_s_a);
% hold on
% plot(reaching_dff_mean_neg_s_a);   


% to indicate the NRAM cell active or not
NRAM_pos_ind_a = [];
NRAM_neg_ind_a = [];
count_pos = 1;
for test_pos = 1:length(NRAM_pos_ind)
    if active_ind(NRAM_pos_ind(test_pos)) == 1
        NRAM_pos_ind_a(count_pos) = NRAM_pos_ind(test_pos);
        count_pos = count_pos+1;
    end
end

count_neg = 1;
for test_neg = 1:length(NRAM_neg_ind)
    if active_ind(NRAM_neg_ind(test_neg)) == 1
        NRAM_neg_ind_a(count_neg) = NRAM_neg_ind(test_neg);
        count_neg = count_neg+1;
    end
end
        
        


%% to calculate mean during whole trial
for num_trial = 1:size(df_Ff_success, 2)
    df_Ff_mean = df_Ff_success(num_trial).df_Ff_mean;
    df_Ff_mean_pos(num_trial) = nanmean(df_Ff_mean(NRAM_pos_ind));
    df_Ff_mean_neg(num_trial) = nanmean(df_Ff_mean(NRAM_neg_ind));
end
data_RB_calcium(session_num).dff_mean_trial.df_Ff_mean_trial_pos = df_Ff_mean_pos';
data_RB_calcium(session_num).dff_mean_trial.df_Ff_mean_trial_neg = df_Ff_mean_neg';

% figure
% plot(df_Ff_mean_pos);
% hold on
% plot(df_Ff_mean_neg); 


%% to calculate mean during the whole session
df_session = data(session_num).df_Ff;
df_session_mean = nanmean(df_session,2);
df_session_pos = df_session_mean(NRAM_pos_ind);
df_session_neg = df_session_mean(NRAM_neg_ind);


%% to calculate heatmap
num_trial = 1;  %% change trial to see a good representative
real_trial = df_Ff_success(num_trial).trail;
active_cell = df_Ff_success(num_trial).active_cell;
reaching_dff = df_Ff_success(num_trial).reaching_dff;
for num_cell = 1:length(active_cell);
    if active_cell(num_cell) == 0
       reaching_dff(num_cell,:) = NaN;
    end
end
dff_max = max(reaching_dff(:));
reaching_dff_norm = reaching_dff/dff_max;
reaching_dff_norm_pos =  reaching_dff_norm(NRAM_pos_ind,:); 
reaching_dff_norm_neg =  reaching_dff_norm(NRAM_neg_ind,:);

% for NRAM_pos SOM
[M I] = max(reaching_dff_norm_pos,[], 2);
max_ind_temp = [I NRAM_pos_ind'];
max_ind = [];
for num_cell = 1:length(M)
    if isnan(M(num_cell))
        max_ind = max_ind;
    else
        max_ind = [max_ind; max_ind_temp(num_cell,:)];
    end
end
% sort the two columns in the same order        
max_ind_x = sortrows(max_ind,1); % based on column #1
cell_peak = max_ind_x(:,2);
dff_heatmap_pos = [];
for num_cell = 1: length(cell_peak)
    dff_heatmap_pos = [dff_heatmap_pos; reaching_dff_norm(cell_peak(num_cell),:)];
end

% colormap('hot')
clims = [0 0.15];
imagesc(dff_heatmap_pos,clims)
colorbar
xticks([])
yticks([])
title_text = [mouse_name ' ' T_Day ' Trial #' num2str(real_trial)];
title(title_text);
            

% for NRAM_neg SOM          
[M I] = max(reaching_dff_norm_neg,[], 2);
max_ind_temp = [I NRAM_neg_ind'];
max_ind = [];
for num_cell = 1:length(M)
    if isnan(M(num_cell))
        max_ind = max_ind;
    else
        max_ind = [max_ind; max_ind_temp(num_cell,:)];
    end
end
% sort the two columns in the same order        
max_ind_x = sortrows(max_ind,1); % based on column #1
cell_peak = max_ind_x(:,2);
dff_heatmap_neg = [];
for num_cell = 1: length(cell_peak)
    dff_heatmap_neg = [dff_heatmap_neg; reaching_dff_norm(cell_peak(num_cell),:)];
end

% colormap('hot')
clims = [0 0.15];
imagesc(dff_heatmap_neg,clims)
colorbar
xticks([])
yticks([])
title_text = [mouse_name ' ' T_Day ' Trial #' num2str(real_trial)];
title(title_text);
            
 

heatmap_pos = reaching_dff_norm_pos(any(~isnan(reaching_dff_norm_pos),2),:);
heatmap_neg = reaching_dff_norm_neg(any(~isnan(reaching_dff_norm_neg),2),:);

% % for NRAM_pos SOM
% [M I] = max(heatmap_pos,[], 2);


%% to calculate the ITI
% to calculate the ITI frame start and end
for n = 1:length(df_Ff_success)  %df_Ff_new is the dff aligned with wavesurfer and calcium data
    temp_frame = df_Ff_success(n).start_Frame;
    if temp_frame < ITI_frame
        df_Ff_success(n).start_ITI_frame = NaN;
        df_Ff_success(n).end_ITI_frame = NaN;
        df_Ff_success(n).dff_ITI = NaN;
        df_Ff_success(n).dff_ITI_mean = NaN;
    else        
        df_Ff_success(n).start_ITI_frame = temp_frame - ITI_frame ;
        df_Ff_success(n).end_ITI_frame = df_Ff_success(n).start_Frame -ITI_frame_end ;
        df_Ff_success(n).dff_ITI = df_Ff_new(:,df_Ff_success(n).start_ITI_frame:df_Ff_success(n).end_ITI_frame);  
        df_Ff_success(n).dff_ITI_mean = nanmean(df_Ff_success(n).dff_ITI,2);
    end
end

% to calculate mean during ITI
for num_trial = 1:size(df_Ff_success, 2)
    if isnan(df_Ff_success(num_trial).dff_ITI_mean)   
        ITI_dff_mean_pos(num_trial) = NaN;
        ITI_dff_mean_neg(num_trial) = NaN;
    else
        ITI_dff_mean = [];
        ITI_dff_mean = df_Ff_success(num_trial).dff_ITI_mean;
        ITI_dff_mean_pos(num_trial) = nanmean(ITI_dff_mean(NRAM_pos_ind_a));
        ITI_dff_mean_neg(num_trial) = nanmean(ITI_dff_mean(NRAM_neg_ind_a)); 
        ITI_dff_mean_total(num_trial) = nanmean(ITI_dff_mean([NRAM_pos_ind_a NRAM_neg_ind_a]));
    end    
end
% reaching_dff_mean_pos_s = sort(reaching_dff_mean_pos);
% reaching_dff_mean_neg_s = sort(reaching_dff_mean_neg);

data_RB_calcium(session_num).ITI_dff_mean.ITI_dff_mean_pos = ITI_dff_mean_pos';
data_RB_calcium(session_num).ITI_dff_mean.ITI_dff_mean_neg = ITI_dff_mean_neg';
data_RB_calcium(session_num).ITI_dff_mean.ITI_dff_mean_total = ITI_dff_mean_total';

%% to calculate failed trial
%to calculate the fail trial
trial_success = [];
for n = 1:length(data_video_success)
    temp = data_video_success(n).trail;
    trial_success = [trial_success temp];
end
trial_success = trial_success;
trial_all = [1:60];
trial_fail = setdiff(trial_all,trial_success);

% set the tone
df_Ff_fail = struct;
for n = 1:length(trial_fail)
    df_Ff_fail(n).trial = trial_fail(n);
    df_Ff_fail(n).tone = data_video(trial_fail(n)).Tone;
    Time_tone_sec = strsplit(df_Ff_fail(n).tone,':');
    Time_tone_start = str2double(Time_tone_sec(1)) * 60 + str2double(Time_tone_sec(2));
    df_Ff_fail(n).tone_frame = round(Time_tone_start*frame_rate);   
    if df_Ff_fail(n).tone_frame ==0
        df_Ff_fail(n).tone_frame =1;
    end
    df_Ff_fail(n).tone_frame_end = round(df_Ff_fail(n).tone_frame + Fail_dura*frame_rate);
    df_Ff_fail(n).dff_fail_cut = df_Ff_new(:,df_Ff_fail(n).tone_frame:df_Ff_fail(n).tone_frame_end);
    mean_fail_cut = nanmean(df_Ff_fail(n).dff_fail_cut,2);
    df_Ff_fail(n).mean_fail_cut = mean_fail_cut;
    mean_fail_cut_pos(n) = nanmean(mean_fail_cut(NRAM_pos_ind));
    mean_fail_cut_neg(n) = nanmean(mean_fail_cut(NRAM_neg_ind));
    
    % for dff before tone (ITI)
    if n>=2
        df_Ff_fail(n).dff_fail_cut_ITI = df_Ff_new(:,df_Ff_fail(n).tone_frame-ITI_frame:df_Ff_fail(n).tone_frame);
        mean_fail_cut_ITI = nanmean(df_Ff_fail(n).dff_fail_cut_ITI,2);
        df_Ff_fail(n).mean_fail_cut_ITI = mean_fail_cut_ITI;
        mean_fail_cut_ITI_pos(n) = nanmean(mean_fail_cut_ITI(NRAM_pos_ind));
        mean_fail_cut_ITI_neg(n) = nanmean(mean_fail_cut_ITI(NRAM_neg_ind));
    end  
%     mean_fail_cut_ITI_pos(n) = mean(mean_fail_cut_ITI(NRAM_pos_ind));
%     mean_fail_cut_ITI_neg(n) = mean(mean_fail_cut_ITI(NRAM_neg_ind));
%         
    
end


        
data_RB_calcium(session_num).fail_dff_mean.mean_fail_cut_pos = mean_fail_cut_pos';
data_RB_calcium(session_num).fail_dff_mean.mean_fail_cut_neg = mean_fail_cut_neg';
data_RB_calcium(session_num).fail_dff_mean.mean_fail_cut_ITI_pos = mean_fail_cut_ITI_pos';
data_RB_calcium(session_num).fail_dff_mean.mean_fail_cut_ITI_neg = mean_fail_cut_ITI_neg';
    

%% to calculate the peak during RB
%to calculate the eligible peak
temp_dff = [];
peak_ind = struct;

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
%         plot_path = 'C:\Users\xyin2\Documents\MATLAB\NPAS4\trace analysis\plot\XY303\';
%         plot(dff_cell);
%         h = hline(threshold_cell,'r');
%         title_text = ['Cell #' num2str(z)];
%         title(title_text);
%         saveas(gcf, [plot_path 'cell #' num2str(z) '.jpg'] );  
%         hold off    
%         close all
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
                        
                        % find three 0 in a row
                        bina_rest = dff_cell_bina(peak_search+1:end);  
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
                        
                        
                        peak_ind(n).trial_num = n;
                        peak_ind(n).frames = nn;                        
                        peak_ind(n).peak(z).cell_num = z;
                        peak_ind(n).peak(z).peak(peak_count).peak_count = peak_count;
                        peak_ind(n).peak(z).peak(peak_count).peak_start = peak_start;
                        peak_ind(n).peak(z).peak(peak_count).peak_end = peak_end;
                        peak_ind(n).peak(z).peak(peak_count).peak_duration = peak_end - peak_start +1;
%                         peak_ind(n).peak(z).peak_count = peak_count;
%                         peak_ind(n).peak(z).peak_start = peak_start;
%                         peak_ind(n).peak(z).peak_end = peak_end;                        

                        dff_event = dff_cell(peak_start:peak_end);
                        peak_ind(n).peak(z).peak(peak_count).peak_max_dff = max(dff_event);
                        
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
                peak_ind(n).trial_num = n;
                peak_ind(n).frames = nn;
                peak_ind(n).peak(z).cell_num = z;                
                peak_ind(n).peak(z).peak(peak_count).peak_count = NaN;
                peak_ind(n).peak(z).peak(peak_count).peak_start = NaN;
                peak_ind(n).peak(z).peak(peak_count).peak_end = NaN;
                peak_ind(n).peak(z).peak(peak_count).peak_duration = NaN;
                peak_ind(n).peak(z).peak(peak_count).peak_max_dff = NaN;
                
                
%                 peak_ind(n).peak(z).peak_count = NaN;
%                 peak_ind(n).peak(z).peak_start = NaN;
%                 peak_ind(n).peak(z).peak_end = NaN;
            end
            
            % to calculate the width and frequency
            cal_table = struct2table(peak_ind(n).peak(z).peak);
            cal_mean = nanmean(cal_table.peak_duration);
            cal_max = nanmax(cal_table.peak_duration); 
            cal_dff_max = nanmax(cal_table.peak_max_dff);
            cal_dff_max_mean = nanmean(cal_table.peak_max_dff);
            peak_ind(n).peak(z).duration = cal_mean;  
            peak_ind(n).peak(z).duration_max = cal_max;  
            peak_ind(n).peak(z).duration_sec = cal_mean/frame_rate;
            peak_ind(n).peak(z).duration_max_sec = cal_max/frame_rate;
            peak_ind(n).peak(z).peak_max_dff = cal_dff_max;
            peak_ind(n).peak(z).peak_max_dff_mean = cal_dff_max_mean;
            
            if isnan(cal_mean)
                peak_ind(n).peak(z).events = NaN;
                peak_ind(n).peak(z).frequency_sec = NaN;
            else               
                peak_ind(n).peak(z).events = length(peak_ind(n).peak(z).peak);
                peak_ind(n).peak(z).frequency_sec = peak_ind(n).peak(z).events/nn*frame_rate;
            end
            
            if ismember(z, NRAM_pos_ind)
                peak_ind(n).peak(z).NRAM_pos = 1;                
            else
                peak_ind(n).peak(z).NRAM_pos = 0;
            end
            
            
            
%               % seperate NRAM_neg vs. NRAM_pos
%                if ismember(z,NRAM_pos)
%                     peak_ind(n).peak_NRAM_pos(NP).cell_num = z;
%                     peak_ind(n).peak_NRAM_pos(NP).peak = peak_ind(n).peak(z).peak;
%                     NP = NP+1;
%                else
%                     peak_ind(n).peak_NRAM_neg(NN).cell_num = z;
%                     peak_ind(n).peak_NRAM_neg(NN).peak = peak_ind(n).peak(z).peak;
%                     NN = NN+1;
%                end

            
            
            
        end
        
%        % seperate NRAM_neg vs. NRAM_pos
% 
%        if ismember(z,NRAM_pos_ind)
%             peak_ind(n).peak_NRAM_pos(NP).cell_num = z;
%             peak_ind(n).peak_NRAM_pos(NP).peak = peak_ind(n).peak(z).peak;
%             NP = NP+1;
%        else
%             peak_ind(n).peak_NRAM_neg(NN).cell_num = z;
%             peak_ind(n).peak_NRAM_neg(NN).peak = peak_ind(n).peak(z).peak;
%             NN = NN+1;
%        end
              
        
    end
end


% plots
duration_total = [];
duration_max_total = [];
frequency_total = [];
peak_max_dff_total = [];
peak_max_dff_mean_total = [];
for xx = 1: length(peak_ind)
    cal_table = [];
    cal_table = struct2table(peak_ind(xx).peak);    
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

data_RB_calcium(session_num).event.peak_ind = peak_ind;
data_RB_calcium(session_num).event.duration_mean_pos = duration_mean_pos ;
data_RB_calcium(session_num).event.duration_mean_neg =duration_mean_neg;
data_RB_calcium(session_num).event.duration_max_pos = duration_max_pos;
data_RB_calcium(session_num).event.duration_max_neg = duration_max_neg;      
data_RB_calcium(session_num).event.frequency_mean_pos = frequency_mean_pos;
data_RB_calcium(session_num).event.frequency_mean_neg = frequency_mean_neg;
data_RB_calcium(session_num).event.peak_max_dff_pos = peak_max_dff_pos;
data_RB_calcium(session_num).event.peak_max_dff_neg = peak_max_dff_neg;
data_RB_calcium(session_num).event.peak_max_dff_pos = peak_max_dff_mean_pos;
data_RB_calcium(session_num).event.peak_max_dff_neg = peak_max_dff_mean_neg;

        
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
            cal_table2 = struct2table(peak_ind(n).peak(z).peak);        
            prob_duration_pos = [prob_duration_pos;cal_table2.peak_duration];
            prob_max_dff_pos = [prob_max_dff_pos;cal_table2.peak_max_dff];
            prob_freq_pos = [prob_freq_pos;peak_ind(n).peak(z).frequency_sec];
        else
            cal_table2 = [];
            cal_table2 = struct2table(peak_ind(n).peak(z).peak);        
            prob_duration_neg = [prob_duration_neg;cal_table2.peak_duration];
            prob_max_dff_neg = [prob_max_dff_neg; cal_table2.peak_max_dff];
            prob_freq_neg = [prob_freq_neg; peak_ind(n).peak(z).frequency_sec];
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
x_pos_prob = sort(prob_duration_pos)/frame_rate;
y_pos_prob = [1:length(prob_duration_pos)]/length(prob_duration_pos);
plot(x_pos_prob, y_pos_prob,'r')

hold on
x_neg_prob = sort(prob_duration_neg)/frame_rate;
y_neg_prob = [1:length(prob_duration_neg)]/length(prob_duration_neg);
plot(x_neg_prob, y_neg_prob,'k')
title_text = ['Event Duration'];
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
title_text = ['Event max dff'];
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
title_text = ['Event frequency'];
title(title_text);
                
 
%% to calculate the low activity cell 
%%(After compare with randomly selected (1,000) NRAM- cells, if it is below 5%, then it is a ‘low activity’ cell)
mean_cell_RB = [];
for num_trial = 1:size(df_Ff_success, 2)
    mean_cell_RB = [mean_cell_RB df_Ff_success(num_trial).reaching_dff_mean];
end
mean_cell_RB_total = nanmean(mean_cell_RB,2);
mean_cell_RB_NRAM_pos = mean_cell_RB_total(NRAM_pos_ind);
mean_cell_RB_NRAM_neg = mean_cell_RB_total(NRAM_neg_ind);
%%%%%%%%
nrand1 = round(length(mean_cell_RB_NRAM_neg)/2);
nrand = 1000;   
%%%%%%%%
low_activity_index = zeros(1,length(mean_cell_RB_NRAM_pos))';
for pos_cell = 1:length(mean_cell_RB_NRAM_pos)
    baseline_pool = [];
    NRAM_neg_cell_num = length(mean_cell_RB_NRAM_neg);
    for rand_x = 1:nrand 
        rand_cell = [];
        rand_cell = mean_cell_RB_NRAM_neg(randi(NRAM_neg_cell_num, nrand1,1));
        difference_pos_neg = nanmean(rand_cell)-mean_cell_RB_NRAM_pos(pos_cell);
        baseline_pool = [baseline_pool difference_pos_neg];
    end
    if prctile(baseline_pool,5)>0
        low_activity_index(pos_cell) =1;
    end
end
per_cell = sum(low_activity_index)/length(mean_cell_RB_NRAM_pos);
disp_M1 = ['low activity cell per. = ' num2str(per_cell)];
disp(disp_M1)


% to locate the NRAM+ cells, to match cells from different session
cell_index_match_session = [];
for cell_check =1:length(low_activity_index)
    cell_index_match_session = [cell_index_match_session find(cell_map_curr==NRAM_pos_ind(cell_check))];
end

% sort cell_index_match_session', but need to relate values in low_activity_index
[Xsorted,I] = sort(cell_index_match_session');
low_activity_index_sorted = low_activity_index(I);
data_output = [Xsorted low_activity_index_sorted];
disp_M2 = ['pie chart data stored in data_output'];
disp(disp_M2)    

disp_M3 = ['major data (mean dff, reliability etc.) stored in data_RB_calcium'];
disp(disp_M3) 

%% control for low activity cells by using the NRAM+ cells 
%%(After compare with randomly selected (1,000) NRAM+ cells, if it is higher 95%, then it is a ‘normal activity’ cell)
nrand2 = round(length(mean_cell_RB_NRAM_pos)/2);
normal_activity_index = zeros(1,length(mean_cell_RB_NRAM_neg))';
for neg_cell = 1:length(mean_cell_RB_NRAM_neg)
    baseline_pool_neg = [];
    NRAM_pos_cell_num = length(mean_cell_RB_NRAM_pos);
    for rand_xx = 1:nrand 
        rand_cell = [];
        rand_cell = mean_cell_RB_NRAM_pos(randi(NRAM_pos_cell_num, nrand2,1));
        difference_pos_neg2 = mean_cell_RB_NRAM_neg(neg_cell)-mean(rand_cell);
        baseline_pool_neg = [baseline_pool_neg difference_pos_neg2];
    end
    if prctile(baseline_pool_neg,95)>0
        normal_activity_index(neg_cell) =1;
    end
end
per_cell = sum(normal_activity_index)/length(mean_cell_RB_NRAM_neg);
disp_M4 = ['normal activity cell per. = ' num2str(per_cell)];
disp(disp_M4)
disp_M5 = ['normal activity cell number = ' num2str(sum(normal_activity_index))];
disp(disp_M5)



NRAM_RB_event_cal_threshold_mean_v2



