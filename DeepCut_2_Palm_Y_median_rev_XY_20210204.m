% read .m files
close all
clear
clc

% mouse name
mouse_name = 'XY311';
% training day
T_Day = 'PR07';
% csv folder
csv_path = 'C:\Users\xyin2\Desktop\XY_Science_2021\csv_data_DLC\';
% deeplabcut_data
DLC_file = [csv_path mouse_name '\' mouse_name '_' T_Day '.mat'];
load(DLC_file)

num_trails = length(data);


% conver string data from all trails into numbers
% calculate median for Head Fixer and Pellet Holder
data_Palm = struct('Head_Fixer_X',[],'Head_Fixer_Y',[],'Palm_X',[],'Palm_Y',[],'Pellet_Holder_X',[],'Pellet_Holder_Y',[], 'Score', [],'HF_Y_median',[],'P_Y_median',[],'PH_Y_median',[]);

for i = 1:num_trails   
    HF_Y = [];
    HF_X = [];
    P_Y = [];   
    P_X = [];
    PH_Y = [];
    PH_X = [];
    Head_F_Y = data(i).head_fixer(:,2);
    Head_F_X = data(i).head_fixer(:,1);
    Palm_Y = data(i).palm(:,2);
    Palm_X = data(i).palm(:,1);
    P_Holder_Y = data(i).pellet_holder(:,2);
    P_Holder_X = data(i).pellet_holder(:,1);
    for j = 1:length(Palm_Y) %%convert string to number
        temp_HF_X = [];
        temp_HF_Y = [];
        temp_P_X = [];
        temp_P_Y = [];
        temp_PH_X = [];
        temp_PH_Y = [];
        temp_HF_X = str2num(Head_F_X{j,1});
        temp_HF_Y = str2num(Head_F_Y{j,1});
        temp_P_X = str2num(Palm_X{j,1});
        temp_P_Y = str2num(Palm_Y{j,1});
        temp_PH_X = str2num(P_Holder_X{j,1});
        temp_PH_Y = str2num(P_Holder_Y{j,1});
        HF_Y(j) = temp_HF_Y;
        HF_X(j) = temp_HF_X;
        P_Y(j) = temp_P_Y;
        P_X(j) = temp_P_X;
        PH_X(j) = temp_PH_X;
        PH_Y(j) = temp_PH_Y;
    end
    data_Palm(i).Head_Fixer_X = HF_X;
    data_Palm(i).Head_Fixer_Y = HF_Y;
    data_Palm(i).Palm_X = P_X;
    data_Palm(i).Palm_Y = P_Y;
    data_Palm(i).Pellet_Holder_X = PH_X;
    data_Palm(i).Pellet_Holder_Y = PH_Y;    
    data_Palm(i).Score = data(i).Score;
    HF_Y_median = median(sort(HF_Y));
    P_Y_median = median(sort(P_Y));
    PH_Y_median = median(sort(PH_Y));
    data_Palm(i).HF_Y_median = HF_Y_median;
    data_Palm(i).P_Y_median = P_Y_median;
    data_Palm(i).PH_Y_median = PH_Y_median;
end

n = 1; % counter of sucssess trails

% data_Palm_sucssess = struct('Trail',[],'LH_X',[],'LH_Y',[],'J_X',[],'J_Y',[], 'Score', [],'JY_median',[]);
data_Palm_sucssess = struct('Trail',[],'Head_Fixer_X',[],'Head_Fixer_Y',[],'Palm_X',[],'Palm_Y',[],'Pellet_Holder_X',[],'Pellet_Holder_Y',[], 'Score', [],'HF_Y_median',[],'P_Y_median',[],'PH_Y_median',[]);

for trail = 1:num_trails
    Score = data_Palm(trail).Score;
    if Score == 11         
        data_Palm_sucssess(n).Trail = trail;
        data_Palm_sucssess(n).Head_Fixer_X = data_Palm(trail).Head_Fixer_X;
        data_Palm_sucssess(n).Head_Fixer_Y = data_Palm(trail).Head_Fixer_Y;
        data_Palm_sucssess(n).Palm_X = data_Palm(trail).Palm_X;
        data_Palm_sucssess(n).Palm_Y = data_Palm(trail).Palm_Y;
        data_Palm_sucssess(n).Pellet_Holder_X = data_Palm(trail).Pellet_Holder_X;
        data_Palm_sucssess(n).Pellet_Holder_Y = data_Palm(trail).Pellet_Holder_Y;  
        data_Palm_sucssess(n).Score = data_Palm(trail).Score;
        data_Palm_sucssess(n).HF_Y_median = data_Palm(trail).HF_Y_median;
        data_Palm_sucssess(n).P_Y_median = data_Palm(trail).P_Y_median;
        data_Palm_sucssess(n).PH_Y_median = data_Palm(trail).PH_Y_median;
        n = n + 1;
    end
end

for i = 1:length(data_Palm_sucssess)
    temp = data_Palm_sucssess(i).Palm_Y;
    data_Palm_sucssess(i).P_Y_max = max(temp);
    data_Palm_sucssess(i).P_Y_min = min(temp);
    data_Palm_sucssess(i).P_Y_sort_median= median(sort(temp));
end

disp ('Wait! Saving the data')
save(save_filename, 'data','mouse_name' ,'T_Day','save_filename','data_Palm','data_Palm_sucssess')
disp ('DONE with Saving')

close all
clear
clc
       
disp ('Done Saving')