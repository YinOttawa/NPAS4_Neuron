%% demo_Irina_part_3
% convert cell_to_index_map into Suite2P_cell_to_index_map
% re-arrange of df/Ff for active ROI's and not active ROI's
% store new data into new_data

disp('running demo_Irina_part_3')
%%

cell_to_index_map = cell_registered_struct.cell_to_index_map;
[n m] = size(cell_to_index_map);
Suite2P_cell_to_index_map = zeros(n,m);

for num_days = 1:m
    temp_Suite2P_cell = (data(num_days).iscell)';
    temp_cell = cell_to_index_map(:,num_days);
    Suite2P_cells = zeros(n,1);
    for i = 1:length(temp_cell)
        temp_S2P_cell = temp_cell(i);
        if temp_S2P_cell == 0
            Suite2P_cells(i) = 0;
        else
            Suite2P_cells(i) = temp_Suite2P_cell(temp_S2P_cell);
        end
    end
    Suite2P_cell_to_index_map(:,num_days) = Suite2P_cells;
end

new_data = struct;
[n, m] = size(cell_to_index_map); % total n - number of ROI's , m - number of Day's
[l, k] = size(data(1).df_Ff);
new_cell_index_map = struct;
temp_sort_cell_to_index_map = cell_to_index_map;
temp_S2P_sort_cell_to_index_map = Suite2P_cell_to_index_map;

for i = 1:m
    temp_sort_cell_to_index_map = sortrows(temp_sort_cell_to_index_map,i);
    temp_S2P_sort_cell_to_index_map = sortrows(temp_S2P_sort_cell_to_index_map,i);
    temp_sort_cell = [];
    temp_sort_S2P_cell = [];
    row = [];
    for j = 1:n
        if temp_sort_cell_to_index_map(j,i) ~= 0
            temp_sort_cell = [temp_sort_cell; temp_sort_cell_to_index_map(j,:)];
            temp_sort_S2P_cell = [temp_sort_S2P_cell; temp_S2P_sort_cell_to_index_map(j,:)];
            row = [row; j];
        end
    end
    new_cell_index_map(i).sort_cells = temp_sort_cell;
    new_cell_index_map(i).S2P_sort_cells = temp_sort_S2P_cell;
    temp_sort_cell_to_index_map = temp_sort_cell_to_index_map(1:row(1)-1,1:m);
    temp_S2P_sort_cell_to_index_map = temp_S2P_sort_cell_to_index_map(1:row(1)-1,1:m);
    [n, m] = size(temp_sort_cell_to_index_map);
end

g = 0;
for i = 1:m
    temp_sort_cell = [];
    temp_sort_cell = new_cell_index_map(i).sort_cells;
    [n m] = size(temp_sort_cell);
    for j = 1:n
        temp_sort_cell_to_index_map = [];
        temp_sort_cell_to_index_map = temp_sort_cell(j,:);
        for l = i:m
            temp_cell = temp_sort_cell_to_index_map(l);
            if temp_cell == 0
                temp_df_Ff = zeros(1,k);
            else
                temp_df_Ff = [];
                temp_df_Ff = data(l).df_Ff(temp_cell,:);
            end
            new_data(l).df_Ff(g+j,:) = temp_df_Ff;
        end
    end
    g = g + j;
end

df_Ff_out_filename = join([df_directory mouse_name '.mat']);
save (df_Ff_out_filename,'data','new_data','new_cell_index_map','cell_to_index_map','Suite2P_cell_to_index_map');
disp('DONE demo_Irina_part_3')   
disp('DONE')   



