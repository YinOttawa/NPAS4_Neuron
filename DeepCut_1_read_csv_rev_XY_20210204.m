% read .csv files and store them into .m files
close all
clear
clc

% mouse name
mouse_name = 'XY311';
% training day
T_Day = 'PR07';
% csv folder
csv_path = 'C:\Users\xyin2\Desktop\XY_Science_2021\csv_data_DLC\';
% location of .csv files
% str1 = 'R:\Simon Chen Laboratory\Lab Members\Xuming\temp\To Irina\XY303_PR04_csv\';
% str1 = 'C:\Users\xyin2\Desktop\XY_Science_2021\csv_data_DLC\XY303\XY303_PR04_csv\';
str1 = [csv_path mouse_name filesep mouse_name '_' T_Day '_csv'];
% location of save . mat files
% str2 = 'C:\Users\xyin2\Desktop\XY_Science_2021\csv_data_DLC\XY303\';
save_path = [csv_path mouse_name '\'];


% select .csv folder and store .csv files names into N array
folder_name = uigetdir(str1);% select folder with .csv files
files = dir(fullfile(folder_name,'*.csv'));
N = natsortfiles({files.name});

% select score data file and store Score data into score_data array
[file_score,path_score] = uigetfile([csv_path mouse_name '\' '*.xlsx']);
file_score_name = [path_score '\' file_score];
file_score_data = readtable(file_score_name);
[m_score n_score] = size(file_score_data);
temp_score_data = table2array(file_score_data);
score_data = temp_score_data(1:m_score,1:n_score);

% number of .csv files
num_files = length(files);
% create struct 
data = struct('head_fixer',[],'palm',[],'pellet_holder',[], 'Score', []);

for t = 1:num_files    
    f_name = char(N(t));
    trail = t;    
    file_name = [folder_name '\' f_name];
    file_data = readtable(file_name);
    [m n] = size(file_data);
    temp_data = table2array(file_data);
    temp_temp_data = temp_data(3:m,2:n);
 
    % Left Hand X & Y direction
    data(trail).head_fixer = temp_temp_data(:,1:2);
 
    % Head Bar X & Y direction
    data(trail).palm = temp_temp_data(:,4:5);
   
    % Joystick X direction
    data(trail).pellet_holder = temp_temp_data(:,7:8);
   
    disp (['Trail #' num2str(trail)])
    data(trail).Score = score_data(t,2);
    
%     prompt = 'What is the score of trail?';
%     % 1 is succsess; 2 is miss; 3 is no attempt
%     x = input(prompt);
%     data(trail).Score = x;
end

disp ('Wait! Saving the data')
save_filename = [save_path mouse_name '_' T_Day '.mat'];
save(save_filename, 'data','mouse_name' ,'T_Day','save_filename')
disp ('DONE with Saving')

% close all
% clear
% clc