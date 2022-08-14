% calculate smooth line trajectory for Sucssess Data_Palm_Sucssess - curv_F
close all
clear
clc

% mouse name
mouse_name = 'XY303';
% training day
T_Day = 'PR04';
% csv folder
csv_path = 'C:\Users\xyin2\Desktop\XY_Science_2021\csv_data_DLC\';
save_path = [csv_path mouse_name '\plot\' 'test\'];
temp_file = [save_path mouse_name '_' T_Day '.mat'];
% deeplabcut_data
DLC_file = [csv_path mouse_name '\' mouse_name '_' T_Day '.mat'];
load(DLC_file)
% load('XY303_PR04.mat')

% pellet holder rich up position
temp = [];
delta_temp = 20;
% % Y value first time passess MAX - 200 and MIN - 100
% % This is start value of movemement   
temp_max_Y = 150;
temp_min_Y = 100;
Palm_Y_xmax = 400;
Palm_Y_xmin = 100;

for i = 1:length(data_Palm_sucssess)
    Palm_Y = [];
    Palm_X = [];
    PH_Y = [];
    PH_Y_median = 0;
    
    Palm_Y = data_Palm_sucssess(i).Palm_Y;% Palm Y
    Palm_X = data_Palm_sucssess(i).Palm_X;% Palm X
    PH_Y = data_Palm_sucssess(i).Pellet_Holder_Y;% Pellet Holder Y
    PH_Y_median = data_Palm_sucssess(i).PH_Y_median;% median of Pellet Holder Y
    num_Palm_Y = length(Palm_Y); 
    time = [1:num_Palm_Y]; % number of frames
    trail_num = data_Palm_sucssess(i).Trail; % Trail number
    pks = findpeaks(Palm_Y);   
    data_Palm_sucssess(i).pks = pks; % peaks of Y
    min_pks = min(pks);
    data_Palm_sucssess(i).min_pks = min_pks; % minimum Peak

    % Pellet holder goes up - PH_Y_frame
    temp_max = PH_Y_median + delta_temp;
    temp_min = PH_Y_median;
    temp = smoothdata(data_Palm_sucssess(i).Pellet_Holder_Y);
    PH_Y_sm = temp;
    for j = 150:length(temp)
        if temp(j) < temp_max
            if temp(j) > temp_min
                temp_frame = j;
                data_Palm_sucssess(i).Pellet_Holder_Y_frame = temp_frame;
                break
            end                 
        end
    end
    
    % Palm Y goes up from min = 100 to max = 150 - PH_Y_frame
    for k = temp_frame:length(Palm_Y)
        if Palm_Y(k) < temp_max_Y
            if Palm_Y(k) > temp_min_Y
                Palm_Y_frame = k;
                data_Palm_sucssess(i).Palm_Y_frame = Palm_Y_frame;
                data_Palm_sucssess(i).f_Palm_Y = Palm_Y(k);
                break
            end                 
        end
    end 
        
    % Plot Palm Y and Peelet Holder Y vs Frames
    fullscreen();
    
    subplot(2,2,1)
    plot (time, Palm_Y,'b');
    hold on        
    plot(time,PH_Y,'k'); 
    hold on
    plot(time,PH_Y_sm,'r'); 
    hold on
    vline(temp_frame,'--r');
    hold on
    vline((Palm_Y_frame - Palm_Y_xmin),'--g');
    hold on
    vline((Palm_Y_frame + Palm_Y_xmax),'--g');
    hold on
    hline((min_pks),'--g');
    hold on
    hline((50),'--r');
    hold on
    xlabel('time (frames)')
    ylabel('Y position in pixel')
    title([mouse_name ' ' T_Day ' Trail ' num2str(trail_num) ' -- Succsess'])
    xmin = 150;
    xmax = 4500;
    xlim([xmin xmax])  
    ylim([0 300])  
    ax = gca;
    ax.YDir = 'reverse';
    
    % Smooth Left Hand Y bottom and top
    Palm_Y_bt =[];
    Palm_X_bt = [];
    
    for m = 1:length(Palm_Y)
        if m > 1
            if Palm_Y(m) > PH_Y_median  % bottom          
                Palm_Y_bt(m) = PH_Y_median; 
                Palm_X_bt(m) = Palm_X(m);  
            elseif Palm_Y(m) < 50  % top 
                if Palm_Y(m-1) > 50
                    temp_Palm_Y_m = Palm_Y(m-1); 
                    temp_Palm_X_m = Palm_X(m-1);               
                    Palm_Y_bt(m) = temp_Palm_Y_m;                
                    Palm_X_bt(m) = temp_Palm_X_m;
                else
                    Palm_Y_bt(m) = temp_Palm_Y_m;                
                    Palm_X_bt(m) = temp_Palm_X_m;         
                end
            else
                Palm_Y_bt(m) = Palm_Y(m); 
                Palm_X_bt(m) = Palm_X(m); 
            end
         end
    end 
    data_Palm_sucssess(i).Palm_Y_bt = Palm_Y_bt;
    data_Palm_sucssess(i).Palm_X_bt = Palm_X_bt;
   
    subplot(2,2,2)
    plot (time, Palm_Y_bt,'b');
    hold on        
%     plot(time,PH_Y,'k'); 
%     hold on
    plot(time,PH_Y_sm,'r'); 
    hold on
    vline(temp_frame,'--r');
    hold on
    vline((Palm_Y_frame - Palm_Y_xmin),'--g');
    hold on
    vline((Palm_Y_frame + Palm_Y_xmax),'--g');
    hold on
    hline((min_pks),'--g');
    hold on
    hline((50),'--r');
    hold on
    xlabel('time (frames)')
    ylabel('Y position in pixel')
    title([mouse_name ' ' T_Day ' Trail ' num2str(trail_num) ' -- Succsess-Smooth'])
    xmin = 150;
    xmax = 4500;
    xlim([xmin xmax])  
    ylim([0 300])  
    ax = gca;
    ax.YDir = 'reverse';
    
    temp_xmin = data_Palm_sucssess(i).Palm_Y_frame - Palm_Y_xmin;
    temp_xmax = data_Palm_sucssess(i).Palm_Y_frame + Palm_Y_xmax;
    if temp_xmax > length(Palm_X)
        temp_xmax = length(Palm_X);
        temp_xmin = length(Palm_X) - (Palm_Y_xmin + Palm_Y_xmax);
    end
    
    P_X = [];
    P_X = data_Palm_sucssess(i).Palm_X(temp_xmin:temp_xmax);
    P_X_bt = [];
    P_X_bt = data_Palm_sucssess(i).Palm_X_bt(temp_xmin:temp_xmax);
    P_Y = [];
    P_Y = data_Palm_sucssess(i).Palm_Y(temp_xmin:temp_xmax);
    P_Y_bt = [];
    P_Y_bt = data_Palm_sucssess(i).Palm_Y_bt(temp_xmin:temp_xmax);
    
    data_Palm_sucssess(i).P_X_bt = P_X_bt;
    data_Palm_sucssess(i).P_Y_bt = P_Y_bt;
    data_P_X_bt = [P_X_bt'];
    data_P_Y_bt = [P_Y_bt'];
    data_P_bt = [data_P_X_bt data_P_Y_bt];
    data_P_bt_sort = sortrows(data_P_bt);
    t = data_P_bt_sort(:,1);
    data_Palm_sucssess(i).t = t;
    y = data_P_bt_sort(:,2);
    data_Palm_sucssess(i).y = y;
    F = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata);
    x0 = [1 1 1 0];
    [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,t,y);
    curv_F = F(x,t);
    data_Palm_sucssess(i).curv_F = curv_F;
    
    PH_X = [];
    PH_X = data_Palm_sucssess(i).Pellet_Holder_X(temp_xmin:temp_xmax);
    PH_Y = [];
    PH_Y = data_Palm_sucssess(i).Pellet_Holder_Y(temp_xmin:temp_xmax);
    
    trail_num = data_Palm_sucssess(i).Trail;
      
    subplot(2,2,3)
    scatter(P_X,P_Y,'b')
    hold on;
    scatter(PH_X,PH_Y,'k')
    hold on;
    title([mouse_name ' ' T_Day ' Trail ' num2str(trail_num) ' -- Succsess'])   
    xlabel('X position in pixel')
    ylabel('Y position in pixel')
    xlim([50 250])  
    ylim([50 250])
    ax = gca;
    ax.YDir = 'reverse';           
    
    subplot(2,2,4)
    scatter(P_X_bt,P_Y_bt,'b')
    hold on;
    scatter(PH_X,PH_Y,'k')
    hold on;
    plot(t,F(x,t),'r')
    hold on;
    title([mouse_name ' ' T_Day ' Trail ' num2str(trail_num) ' -- Succsess-Smooth'])   
    xlabel('X position in pixel')
    ylabel('Y position in pixel')
    xlim([50 250])  
    ylim([50 250])
    ax = gca;
    ax.YDir = 'reverse';
%     str = ['R:\Simon Chen Laboratory\Lab Members\Xuming\temp\To Irina\DeepCut_3_plot\'];
    saveas(gcf, [save_path 'Success Trail_' mouse_name '_TDay_' T_Day ' Trail ' num2str(trail_num) '.jpg'] ); 
    hold off    
    close all
end

% save(filename, 'data','mouse_name','T_Day','filename','data_Palm','data_Palm_sucssess');
save(temp_file, 'data','mouse_name','T_Day','save_filename','data_Palm','data_Palm_sucssess');