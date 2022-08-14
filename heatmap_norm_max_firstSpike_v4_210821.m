
ITI_sec = 10;
ITI_frame = round(ITI_sec*frame_rate);
video_framerate = 150;

root_path = 'C:\Users\xyin2\Documents\MATLAB\NPAS4\trace analysis\';
save_path_heatmap = [root_path mouse_name '\'];

DLC_file = 'F:\Suite2P_analysis\XY305\XY305_PR04\XY305_PR04_DLC';


for trial_num = 1:length(df_Ff_success);% 4,5 good
    real_trial = df_Ff_success(trial_num).trail;

    % sort the cells by the begining of first event (peak), event_start_index_new
    event_start = [];
    for n = 1:length(data_RB_calcium(session_num).event.peak_ind(trial_num).peak)    
        event_start(n) = data_RB_calcium(session_num).event.peak_ind(trial_num).peak(n).peak.peak_start;
    end
    event_start_order = [1:1:length(event_start)];
    event_start_index = [event_start_order' event_start'];
    event_start_new = sortrows(event_start_index,2); 
    event_start_index_new = event_start_new(:,1);

    % locate the RB from the all dff, then include the ITI
    temp_reaching_dff = df_Ff_success(trial_num).reaching_dff;
    temp_reaching_cell = temp_reaching_dff(n,:);
    temp_dff_cell = df_Ff_new(n,:);
    [Lia, Locb] = ismember(temp_dff_cell,temp_reaching_cell);
    RB_start = find(Locb,1);
    RB_end = RB_start + size(temp_reaching_dff,2);
    ITI_start = RB_start - ITI_frame +1;
    temp_dff = df_Ff_new(:,[ITI_start:RB_end]);  
   
    
    dff_heatmap_RB = [];
    % temp_dff = df_Ff_success(trial_num).reaching_dff;
    for n =1:length(event_start_index_new);  
        dff_heatmap_RB(n,:) = temp_dff(event_start_index_new(n),:);
    end
    % colormap('hot')
    % clims = [nanmean(df_baseline_threshold) 2];
    clims = [0 2];
    imagesc(dff_heatmap_RB,clims)
    colorbar
    xticks([0:50:size(dff_heatmap_RB,2)])
    yticks([0:10:size(dff_heatmap_RB,1)])
    % title_text = [mouse_name ' ' T_Day ' Trial #' num2str(real_trial) ' NRAM-Neg'];
    title_text = ['Trial #' num2str(real_trial) ' All SOM'];
    title(title_text);

    hold on
    h = vline(RB_start-ITI_start,'r--','RB start');
    
    % locate the 'pellet to mouth' frame
    load(DLC_file);
    pellet_to_mouth_frame = size(temp_dff,2)-round((data_DLC_ind(trial_num).reaching_end - ...
        data_DLC_ind(trial_num).palm_mouth)/video_framerate*frame_rate);
    
    hold on
    h = vline(pellet_to_mouth_frame,'r--','pellet to mouth');
    
    % locate the 'tone'
    reaching_end_frame = round(data_DLC_ind(trial_num).reaching_end/video_framerate*frame_rate);
    tone_frame = size(temp_dff,2)- reaching_end_frame -round(1*frame_rate); % 1 second tone
    
    if tone_frame > 0
        h = vline(tone_frame,'r--','tone');
    end       

    
    saveas(gcf, [save_path_heatmap mouse_name ' Trial #' num2str(real_trial) '.jpeg'] );  
    hold off    
    close all
end

