close all
clear
clc


% mouse name
mouse_name = 'XY303';
% training day
T_Day = 'PR04';
csv_path = 'C:\Users\xyin2\Desktop\XY_Science_2021\csv_data_DLC\';

save_path = [csv_path mouse_name '\plot\' ];

% deeplabcut_data
DLC_file = [csv_path mouse_name '\' mouse_name '_' T_Day '.mat'];

load(DLC_file)

num_trails = length(data_Palm_sucssess);
figure('Name','Sucssess Trails');

for i = 1:num_trails
    Palm_Y = data_Palm_sucssess(i).Palm_Y;
    num_Palm_Y = length(Palm_Y);
    
%     HF_Y = data_Palm_sucssess(i).Head_Fixer_Y;
    HF_Y_median = data_Palm_sucssess(i).HF_Y_median;
    HF_Y = HF_Y_median * ones(1,num_Palm_Y,'uint16');
    
%     PH_Y = data_Palm_sucssess(i).Pellet_Holder_Y;
    PH_Y_median = data_Palm_sucssess(i).PH_Y_median;
    PH_Y = PH_Y_median * ones(1,num_Palm_Y,'uint16');
   
    Trail = data_Palm_sucssess(i).Trail;
   
    time = [1:num_Palm_Y]; 
    fullscreen();  %%figure fullscreen
    subplot(3,6,i)
    plot (time, Palm_Y,'k');
    hold on        
    plot(time,HF_Y,'g');     
    hold on        
    plot(time,PH_Y,'r');  
    hold on
    xlabel('time (frames)')
    ylabel('Y position in pixel')
    title([mouse_name ' ' T_Day ' Trail ' num2str(Trail)])
    xmin = 150;
    xmax = 4500;
    xlim([xmin xmax])  
    ax = gca;
    ax.YDir = 'reverse';
end


saveas(gcf, [save_path 'DeepCut_2_Plot_Success Trail_' mouse_name '_TDay_' T_Day '_all_trails_' num2str(xmin) '_' num2str(xmax) '.jpg'] ); 
saveas( gcf, [save_path 'DeepCut_2_Plot_Success Trail_' mouse_name '_TDay_' T_Day '_all_trails_' num2str(xmin) '_' num2str(xmax) '.eps'] ,'epsc');
        
