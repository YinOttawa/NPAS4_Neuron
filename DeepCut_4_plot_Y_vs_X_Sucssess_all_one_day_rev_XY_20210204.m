% plot Y vs X Success Trail trajectory
close all
clear
clc

load('XY303_PR04.mat')

fullscreen();
for i = 1:length(data_Palm_sucssess)
    F = data_Palm_sucssess(i).curv_F;
    t = data_Palm_sucssess(i).t;
    plot(t,F,'k')
    hold on;
    title([mouse_name ' ' T_Day])   
    xlabel('X position in pixel')
    ylabel('Y position in pixel')
    xlim([50 250])  
    ylim([50 250])
    ax = gca;
    ax.YDir = 'reverse';    
end
str = ['R:\Simon Chen Laboratory\Lab Members\Xuming\temp\To Irina\'];
saveas(gcf, [str 'DeepCut_4_PLOT_Y_ vs_X_Success Trail_' mouse_name '_TDay_.jpg'] ); 
hold off    
% close all
