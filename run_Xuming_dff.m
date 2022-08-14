%% Please change Mouse Name and Directory:
% Mouse Name
mouse_name = 'XY312';% Xuming example
% mouse_name = 'CL198';% Candice example

% Defining the directory
% directory = 'R:\Simon Chen Laboratory\Codes\Irina\CellReg-master\test_Irina';
directory = 'C:\Users\xyin2\Documents\MATLAB\Suite2P\dff_extraction\Result_output';

output_dir = 'F:\Suite2P_output';

%% The results_directory for CellReg SampleData:
results_directory = [directory '\SampleData\'];

if exist(results_directory,'dir')~=7
    mkdir(results_directory);
end

%% The df_directory for Suite2P Data:
df_directory = [directory '\df_Ff_Data\'];

if exist(df_directory,'dir')~=7
    mkdir(df_directory);
end

%% Loading the Suite2P data:
% Suite2P_file_path = 'E:\S2P\CellReg-master\test_Irina\Sute2PData\';
% Suite2P_file_names = {'F_XY212_2020-01-15_plane1_proc.mat' ,...
%             'F_XY212_2020-01-16_plane1_proc.mat' ,...
%             'F_XY212_2020-01-18_plane1_proc.mat'};

[Suite2P_file_names,Suite2P_file_path] = uigetfile([output_dir '\' '*.mat'],...
   'Select at least TWO or More Suite2P DATA Files', ...
   'MultiSelect', 'on');

if isequal(Suite2P_file_names,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(Suite2P_file_path,Suite2P_file_names{1,1})])
end

%% demo_Irina_part_1 creates images of each Sute2PData and saves as data_out
% into SampleData folder and into variable file_names:
% file_names = {'E:\S2P\CellReg-master\test_Irina\SampleData\XY212_day1_CellReg.mat' ,...
%             'E:\S2P\CellReg-master\test_Irina\SampleData\XY212_day2_CellReg.mat' ,...
%             'E:\S2P\CellReg-master\test_Irina\SampleData\XY212_day3_CellReg.mat'};
% 
% df_Ff data saves into df_Ff_Data folder
% use AP_baselineEstimation.m to calucalte df_Ff and baseline

demo_Irina_part_1_1

%% demo_Irina_part_2 aligning all sessions ROI's and creates cell_to_index_map
%
% modified vertion of CellReg demo
% functions used by CellReg:
% adjust_FOV_size.m
% align_images.m
% check_if_in_overlapping_FOV.m
% choose_best_model.m
% cluster_cells.m
% compute_centroid_distances_model.m
% compute_centroid_locations.m
% compute_centroid_projections.m
% compute_data_distribution.m
% compute_footprint_projections.m
% ompute_p_same.m
% compute_scores.m
% display_progress_bat.m
% estimate_number_of_bins.m
% estimate_registration_accuracy.m
% evaluate_data_quality.m
% gaussfit.m
% initial_registration_centroid_distances.m
% interpolate_pixel_value.m
% load_multiple_sessions.m
% normalized _spatial_footprints.m
% plot_alignment_results.m
% plot_all_registered_projections.m
% plot_all_sessions_projects.m
% plot_cell_scores.m
% plot_estimated_registration_accuracy.m
% plot_initial_registration.m
% plot_models.m
% plot_x_y_displacements.m
% translate_projections.m
%
% Input file_names stored by part1 into SampleData folder and into variable
% file_names

results_directory = [directory '\SampleResults\'];

if exist(results_directory,'dir')~=7
    mkdir(results_directory);
end

demo_XY_part_2

%% demo_Irina_part_3 

demo_XY_part_3

%% THE END
clear all;
close all;
clc;
%%
        