close all
clear
clc

% mouse name
mouse_name = 'XY331';
% training day
T_Day = 'PR05';
% csv folder
csv_path = 'C:\Users\xyin2\Desktop\XY_Science_2021\csv_data_DLC\';
save_path = [csv_path mouse_name '\plot\' 'reaching_bout\' T_Day '\'];
if exist(save_path,'dir')~=7  % 7: name is a folder.
    mkdir(save_path);
end
save_file = [save_path mouse_name '_' T_Day '_DLC''.mat'];
% temp_file = [save_path mouse_name '_' T_Day '.mat'];
% deeplabcut_data
DLC_file = [csv_path mouse_name '\' mouse_name '_' T_Day '.mat'];
load(DLC_file)


% xmin = 150;
% xmax = 4500;

threshold = 120; %% threshold for the mouth height
threshold_PH = 200; %% threshold for the pellet holder height
threshold_duration = 200; %% only mouth time higher than threshold over 200 frmaes would be counted, XY305 trial 31 reaching late 
subtraction_reverse = 352;% video size: 260*352
duration_end = 300; %duration after pellet was sent to mouth 

data_DLC_ind = [];
data_DLC_ind.trial = [];
data_DLC_ind.palm_holder = [];
data_DLC_ind.PH_ready = [];
data_DLC_ind.palm_mouth =  [];
data_DLC_ind.reaching_end = [];
    
for i = 1:length(data_Palm_sucssess)
    Palm_Y = [];
    Palm_Y = data_Palm_sucssess(i).Palm_Y;
    num_Palm_Y = length(Palm_Y); 
    trial_num = data_Palm_sucssess(i).Trail; % Trail number


    %% make binary and find the time when the palm sends the pellet to the mouth
    binary.Palm_Y = zeros(num_Palm_Y,1);
    for n = 1:num_Palm_Y
        if Palm_Y(n) > threshold
            binary.Palm_Y(n) = 0;
        else
            binary.Palm_Y(n) = 1;  %% because the original data is reversed, smaller data is the ones we want
        end
    end
    % when the duration is over 200 frames (1s), the mouth time will be valid
    mouth_pellet_candidate = [];
    for m = 1:num_Palm_Y
        if  binary.Palm_Y(m) == 1;
            temp_duration = [];
            temp_duration = binary.Palm_Y(m:m+threshold_duration);
            if sum(temp_duration) > threshold_duration
                mouth_pellet_candidate = m;
                break
            end
        end
    end
    % exclude the noise data (before a valid mouth sending)
    temp_mouth = [];
    for p = 1: num_Palm_Y;
        if p < mouth_pellet_candidate
            temp_mouth(p) = subtraction_reverse;%%make it excluded from the rest data
        else
            temp_mouth(p) = Palm_Y(p);
        end
    end

    % reverse the data by (352-x) and find the peak (a alternative way to find the minimun)
    temp_mouth = subtraction_reverse-temp_mouth;
    [pks,pks_ind] = findpeaks(temp_mouth);
    mouth_start = pks_ind(1); % the first peak as the first frame to the mouth 
    
    %% find the reaching bout start
    palm_holder = [];
    Pellet_Holder_Y = data_Palm_sucssess(i).Pellet_Holder_Y;
    PH_median = median(Pellet_Holder_Y);
    PH_median_rev = subtraction_reverse - PH_median;
    temp_reaching = [];
    temp_reaching_rev = [];
    temp_reaching = Palm_Y(1:mouth_start);
    temp_reaching_rev = subtraction_reverse - temp_reaching;
    [r_pks,r_pks_ind] = findpeaks(temp_reaching_rev);
    for q = 1:length(r_pks)
        if r_pks(q)> PH_median_rev
            palm_holder = r_pks_ind(q);
            break
        end
    end 
    
    % if palm_holder = [], which means no peak before the mouth_time, then
    % choose the time point passing the 'PH_median_rev'
    if isempty(palm_holder)
        temp_reaching_rev_temp = temp_reaching_rev;
        temp_reaching_rev_temp(temp_reaching_rev_temp<PH_median_rev)=0;
        [temp_find, temp_find_ind] = find(temp_reaching_rev_temp);
        palm_holder = temp_find_ind(1); % the first value passing the 'PH_median_rev'
    end         
    
    
    
    %% find pellet hold reach its height
%     plot (time, Pellet_Holder_Y,'r');
    binary.PH = zeros(num_Palm_Y,1);
    for r = 1:num_Palm_Y
        if Pellet_Holder_Y(r) > threshold_PH
            binary.PH(r) = 0;
        else
            binary.PH(r) = 1;
        end
    end
    % when the duration is over 200 frames (~1s), the pellet time will be valid
    PH_ready = [];
    for s = 1: num_Palm_Y
        if binary.PH(s) == 1
            temp_duration = [];
            temp_duration = binary.PH(s:s+threshold_duration);
            if sum(temp_duration) > threshold_duration
                PH_ready = s;
                break
            end
        end        
    end
%     vline(PH_ready,'--k'); 
     
%     for k = mouth_start:-1:1;
%         if Palm_Y(k) >= Palm_Y(k-1)
%             reaching_start = k;
%             break
%         end
%     end
    %% make database
    data_DLC_ind(i).trial = data_Palm_sucssess(i).Trail;
    data_DLC_ind(i).palm_holder = palm_holder ;
    data_DLC_ind(i).PH_ready = PH_ready;
    data_DLC_ind(i).palm_mouth =  mouth_start;
    reaching_end = mouth_start + duration_end;
    data_DLC_ind(i).reaching_end = reaching_end;
    
    

    %% plot figure
    figure
    time = [1:num_Palm_Y]; % number of frames
    Palm_Y_rev = [];
    Palm_Y_rev = 352 - Palm_Y;  % video size: 260*352
    plot (time, Palm_Y_rev,'b');
    % ax = gca;
    % ax.YDir = 'reverse';
    % xlim([xmin xmax]);
    ylim([0 300]);
    vline(mouth_start,'--r');
    vline(palm_holder,'--k');
    vline(reaching_end,'--k');
    title([mouse_name ' ' T_Day ' Trial ' num2str(trial_num) ' -- Succsess'])
    saveas(gcf, [save_path 'Success Trail_' mouse_name '_TDay_' T_Day ' Trial ' num2str(trial_num) '_test.jpg'] ); 
    saveas(gcf, [save_path 'Success Trail_' mouse_name '_TDay_' T_Day ' Trial ' num2str(trial_num) '.pdf'] ); 
    hold off    
    close all
    

    
end

save(save_file, 'mouse_name','T_Day','data_DLC_ind');





