%% demo_Irina_part_1
% convert Suite2P data into CellReg.mat data
% calculation of df/Ff for active ROI's

disp('running demo_Irina_part_1')
number_of_sessions = size(Suite2P_file_names,2);
data = struct;
for n = 1:number_of_sessions
    this_file_name = fullfile(Suite2P_file_path,Suite2P_file_names{1,n});
    temp_data = load(this_file_name);
    [k m] = size (temp_data.dat.stat); 
    
    % iscell
    iscell = [];
    
    Ff = [];
    Fcell = [];
    FcellNeu = [];
    for i = 1:m
        if temp_data.dat.stat(i).iscell == 1
            iscell = [iscell i];  
            temp_Ff = temp_data.dat.Ff(:,i)';
            temp_Fcell = temp_data.dat.Fcell{1,1}(i,:);
            temp_FcellNeu = temp_data.dat.FcellNeu{1,1}(i,:);
            Ff = [Ff; temp_Ff];
            Fcell = [Fcell; temp_Fcell];
            FcellNeu = [FcellNeu; temp_FcellNeu];  
        end
    end
    data(n).iscell = iscell;
    data(n).Ff = Ff;
    data(n).Fcell = Fcell;
    data(n).FcellNeu = FcellNeu;
    
    % dF/F Calculation
%     framerate = 1;
    framerate = 30; %% XY_modi
%     roi_trace_long = [];
    Ff(isnan(Ff))=0;
    roi_trace_long = Ff;      
    roi_trace_df = AP_baselineEstimation(roi_trace_long,framerate);  
    data(n).df_Ff = roi_trace_df;
    
    Lx = temp_data.dat.cl.Lx;
    Ly = temp_data.dat.cl.Ly;
    data_out = zeros(length(iscell),Ly,Lx);

    for k = 1:length(iscell)
        out = zeros(Ly,Lx);
        ncell = iscell(k);
        npix = temp_data.dat.stat(ncell).npix;
        n_ypix = temp_data.dat.stat(ncell).ypix;
        n_xpix = temp_data.dat.stat(ncell).xpix;
        n_ipix = temp_data.dat.stat(ncell).ipix;
        for i = 1:npix
            out(n_ypix(i),n_xpix(i)) = n_ipix(i);
        end
        data_out(k,:,:) = out;
    end
    data_out_filename = [results_directory mouse_name '_day' num2str(n) '_CellReg.mat'];
    save (data_out_filename,'data_out');  
    disp(['DONE DAY' num2str(n)])  
    file_names(n) = {data_out_filename};
end

df_Ff_out_filename = join([df_directory mouse_name '.mat']);
save (df_Ff_out_filename,'data');
disp('DONE demo_Irina_part_1')   